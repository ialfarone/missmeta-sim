# I simulate here a MNAR scenario for a mvmeta analysis of CR and SR.

# Cond 1: Missingness -> 5%, 10%, 15%, 20%
# Cond 2: Delta adjustment CR -> - 3, - 1.5, 0, 1.5, 3
# Cond 2bis: Delta adjustment SR -> - 3, - 1.5, 0, 1.5, 3
# Cond 3: Distribution -> uniform, normal(0), normal(3), normal(-3), mvt(3,3)

# I have then 4*5*5*5 = 500 situations
# I simulate one dataset, I perform all the analyses and then I replicate them 10^4 times

# I want to extract a list in which I have a dataset for each condition
# The goal is to look at the distribution of the CR and SR

library(mvtnorm)
library(mixmeta)
library(systemfit)
library(dplyr)
library(ggplot2)
library(tmvtnorm)

S = 50
N = 100

Mu = c(0, 0)
Tau = c(1, 1)
rho = 0.7
Sigma = diag(Tau) %*% matrix(c(1, rho, rho, 1), nrow = 2) %*% diag(Tau)
RTher = rmvnorm(S, Mu, Sigma)

b0 = rnorm(S, mean = 20, sd = 2)
b1 = rnorm(S, mean = 0.5, sd = 0.2)
b2 = rnorm(S, mean = 1.5, sd = 0.3)
b3 = 3

data = vector(mode = "list", length = S)

for (i in 1:S) {
  minA = runif(1, min = 18, max = 65)
  maxA = runif(1, min = minA + 5, max = 90)
  Age = runif(N, min = minA, max = maxA)
  
  pFem = 0.45
  Sex = rbinom(N, size = 1, prob = pFem)
  Ther = rbinom(N, size = 1, prob = 0.5)
  
  CR = rnorm(N,
             mean = b0 + b1 * Age + b2 * Sex + (b3 + RTher[i, 1]) * Ther,
             sd = 5)
  SR = rnorm(N,
             mean = b0 + b1 * Age + b2 * Sex + (b3 + RTher[i, 2]) * Ther,
             sd = 9)
  
  data[[i]] = data.frame(
    Study = i,
    Age = Age,
    Sex = factor(Sex, levels = 0:1, labels = c("M", "F")),
    Ther = factor(Ther, levels = 0:1, labels = c("New", "Std")),
    CR = CR,
    SR = SR
  )
}

d = do.call(rbind, data)

dat = vector(mode = "list", length = S)

for (s in 1:S) {
  Sn = d[d$Study == s, ]
  
  m1 = CR ~ Age + Sex + Ther
  m2 = SR ~ Age + Sex + Ther
  
  fitsur = systemfit(list(CR = m1, SR = m2), "SUR", data = Sn)
  sum = summary(fitsur)
  
  dat[[s]] = data.frame(
    Study = s,
    EstCR = sum$coefficients[4, 1],
    EstSR = sum$coefficients[8, 1],
    SECR = sum$coefficients[4, 2],
    SESR = sum$coefficients[8, 2],
    Cor.ws = sum$residCor["CR", "SR"]
  )
}

dat = do.call(rbind, dat)
head(dat)

## Multivariate meta-analysis ####

theta = cbind(dat$EstCR, dat$EstSR)
cor2cov = function(sd1, sd2, rho) {
  sd1 * sd2 * rho
}
Sigma = cbind(dat$SECR^2,
              cor2cov(dat$SECR, dat$SESR, dat$Cor.ws),
              dat$SESR^2)

mv.c = mixmeta(theta, Sigma, method = "reml")
summary(mv.c)

# Generate missing ####

## Missing Not At Random ####
# --- Directional Exclusive MNAR ---
target <- 0.2
betaCR <- 2     # CR missing when low
betaSR <- -1.5  # SR missing when high

invlogit <- function(x) plogis(x)

beta0_cr <- uniroot(function(b0) mean(invlogit(b0 + betaCR * dat$EstCR)) - (1 - target / 2), c(-20, 20))$root
prob_cr <- invlogit(beta0_cr + betaCR * dat$EstCR)

beta0_sr <- uniroot(function(b0) mean(invlogit(b0 + betaSR * dat$EstSR)) - (1 - target / 2), c(-20, 20))$root
prob_sr <- invlogit(beta0_sr + betaSR * dat$EstSR)

M_cr <- rbinom(nrow(dat), 1, 1 - prob_cr)
M_sr <- rbinom(nrow(dat), 1, 1 - prob_sr)

# Conflict resolution: keep the more extreme one
conflict <- which(M_cr == 1 & M_sr == 1)
meanCR <- mean(dat$EstCR)
meanSR <- mean(dat$EstSR)

for (i in conflict) {
  dist_cr <- abs(dat$EstCR[i] - meanCR)
  dist_sr <- abs(dat$EstSR[i] - meanSR)
  if (dist_cr > dist_sr) {
    M_sr[i] <- 0  # CR is more deviant: keep SR, drop CR
  } else {
    M_cr[i] <- 0  # SR is more deviant: keep CR, drop SR
  }
}

dmnar <- dat
dmnar[M_cr == 1, c("EstCR", "SECR")] <- NA
dmnar[M_sr == 1, c("EstSR", "SESR")] <- NA
dmnar$Cor.ws[is.na(dmnar$EstCR) | is.na(dmnar$EstSR)] <- NA
# Assuming `dat` is your study-level data frame with effect estimates

mnar_result <- genmnar(dat, target = 0.2, betaCR = 2, betaSR = -1.5)

# Extract the MNAR dataset
dmnar <- mnar_result$data

# Print diagnostic info
print(mnar_result$table)
cat("Overall missing rate:", round(100 * mnar_result$missing_rate, 2), "%\n")


table(
  CR_missing = is.na(dmnar$EstCR),
  SR_missing = is.na(dmnar$EstSR)
)

ggplot(dmnar, aes(x = dat$EstSR, fill = is.na(EstSR))) +
  geom_histogram(binwidth = 2, position = "stack") +
  labs(title = "Missingness in CR as a function of CR", x = "EstCR", y = "Count") +
  scale_fill_manual(values = c("#5B5F97", "#EE6C4D"), name = "CR Missing") +
  theme_minimal()

## Multivariate meta-analysis #### Check

theta.m = cbind(dmnar$EstCR, dmnar$EstSR)
Sigma.m = cbind(dmnar$SECR^2,
                cor2cov(dmnar$SECR, dmnar$SESR, dmnar$Cor.ws),
                dmnar$SESR^2)

mv.m = mixmeta(theta.m, Sigma.m, method = "reml")
summary(mv.m)


genimp.mnar = function(df, distribution = c("uniform", "normal", "tmvn"),
                       iter, minCR, maxCR, minSR, maxSR, meanCR, meanSR, sdCR, sdSR,
                       meantmv, sigmatmv, lower, upper, imprho, delta, scaleSE) {
  
  distribution = match.arg(distribution)
  results = list()
  
  cor2cov = function(sd1, sd2, rho) {
    sd1 * sd2 * rho
  }
  
  for (d in 1:length(delta)) {
    deltadj = delta[d]  # renamed from delta_value
    
    for (i in 1:iter) {
      dfi = df 
      
      NmissCR = sum(is.na(dfi$EstCR))
      NmissSR = sum(is.na(dfi$EstSR))
      
      if (distribution == "uniform") {
        
        impCR = runif(NmissCR, min = minCR, max = maxCR)
        impSR = runif(NmissSR, min = minSR, max = maxSR)
        
      } else if (distribution == "normal") {
        
        impCR = rnorm(NmissCR, mean = meanCR + deltadj, sd = sdCR)
        impSR = rnorm(NmissSR, mean = meanSR + deltadj, sd = sdSR)
        
      } else if (distribution == "tmvn") {
        
        imputed = rtmvnorm(
          n = sum(is.na(dfi$Cor.ws)),
          mean = meantmv + c(deltadj, deltadj),
          sigma = sigmatmv,
          lower = lower,
          upper = upper
        )
        
        impCR = imputed[, 1]
        impSR = imputed[, 2]
      }
      
      dfi$EstCR[is.na(dfi$EstCR)] = sample(impCR, NmissCR, replace = FALSE)
      dfi$EstSR[is.na(dfi$EstSR)] = sample(impSR, NmissSR, replace = FALSE)
      dfi$Cor.ws[is.na(dfi$Cor.ws)] = imprho
      
      impSECR = sample(dfi$SECR[!is.na(dfi$SECR)], NmissCR, replace = TRUE) * scaleSE
      impSESR = sample(dfi$SESR[!is.na(dfi$SESR)], NmissSR, replace = TRUE) * scaleSE
      
      dfi$SECR[is.na(dfi$SECR)] = impSECR 
      dfi$SESR[is.na(dfi$SESR)] = impSESR 
      
      theta = cbind(dfi$EstCR, dfi$EstSR)
      Sigma = cbind(dfi$SECR^2,
                    cor2cov(dfi$SECR, dfi$SESR, dfi$Cor.ws),
                    dfi$SESR^2)
      
      mv = mixmeta(theta, Sigma, method = "reml")
      ci = confint(mv)
      
      results[[length(results) + 1]] = data.frame(
        deltadj = deltadj,
        iter = i,
        eff1 = mv$coefficients[1],
        eff2 = mv$coefficients[2],
        se1 = sqrt(mv$vcov[1, 1]),
        se2 = sqrt(mv$vcov[2, 2]),
        cov = mv$vcov[1, 2],
        ci.lb1 = ci[1, 1],
        ci.ub1 = ci[1, 2],
        ci.lb2 = ci[2, 1],
        ci.ub2 = ci[2, 2]
      )
    }
  }
  
  do.call(rbind, results)
}

resnorm = genimp.mnar(
  df = dmnar,
  distribution = "normal",
  iter = 10,
  meanCR = 3, 
  meanSR = 3,
  sdCR = 5,
  sdSR = 5,
  imprho = 0.7,
  delta = seq(-7, 7, by = 2),
  scaleSE = 1.5
)

resnorm

##################################
### Function to summarize data####
##################################

sum.meth = function(res, true1, true2, method_name) {
  m = nrow(res)
  
  Q_mat = cbind(res$eff1, res$eff2)
  Q_bar = colMeans(Q_mat)                     
  
  B = cov(Q_mat)                              
  
  U_list = lapply(1:m, function(i) {
    matrix(c(res$se1[i]^2, res$cov[i], 
             res$cov[i], res$se2[i]^2), 
           nrow = 2)
  })
  
  U_bar = Reduce("+", U_list) / m             
  
  Tmat = U_bar + (1 + 1/m) * B # total pooled variance (now I incorporate also the covariance)
  
  se = sqrt(diag(Tmat))
  
  df = (m - 1) * (1 + diag(U_bar) / ((1 + 1/m) * diag(B)))^2
  
  ci1 = Q_bar[1] + c(-1, 1) * qt(0.975, df[1]) * se[1]
  ci2 = Q_bar[2] + c(-1, 1) * qt(0.975, df[2]) * se[2]
  
  bias1 = Q_bar[1] - true1
  bias2 = Q_bar[2] - true2
  
  cover1 = as.numeric(ci1[1] <= true1 && true1 <= ci1[2])
  cover2 = as.numeric(ci2[1] <= true2 && true2 <= ci2[2])
  
  return(data.frame(
    method = method_name,
    est_CR = Q_bar[1],
    est_SR = Q_bar[2],
    bias_CR = bias1,
    bias_SR = bias2,
    cover_CR = cover1,
    cover_SR = cover2,
    pci_lb_CR = ci1[1],
    pci_ub_CR = ci1[2],
    pci_lb_SR = ci2[1],
    pci_ub_SR = ci2[2]
  ))
}

true1 <- 3
true2 <- 3

sumbydelta <- lapply(
  split(resnorm, resnorm$delta),
  function(df) 
    sum.meth(df, true1 = true1, true2 = true2, 
             method_name = paste0("delta_", unique(df$delta)))
)

resbydelta = do.call(rbind, sumbydelta)
resbydelta



ggplot(resbydelta, aes(x = seq(-7, 7, by = 2))) +
  geom_line(aes(y = est_CR), color = "#1f77b4") +
  geom_ribbon(aes(ymin = pci_lb_CR, ymax = pci_ub_CR), alpha = 0.2, fill = "#1f77b4") +
  labs(x = expression(delta),
       y = "Pooled Estimate (CR)") +
  theme_minimal()

ggplot(resbydelta, aes(x = seq(-7, 7, by = 2), y = bias_CR)) +
  geom_line(color = "#d62728") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = expression(delta),
       y = "Bias (Estimate - True Value)") +
  theme_minimal()

ggplot(resbydelta, aes(x = seq(-7, 7, by = 2), y = cover_CR)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen") +
  labs(x = expression(delta),
       y = "Coverage") +
  theme_minimal()


deltas <- seq(-7, 7, by = 2)  # make sure this is defined for x-axis

ggplot(resbydelta, aes(x = deltas)) +
  # CR line and ribbon
  geom_line(aes(y = est_CR), color = "#1f77b4", size = 1) +
  #geom_ribbon(aes(ymin = pci_lb_CR, ymax = pci_ub_CR), alpha = 0.2, fill = "#1f77b4") +

    geom_line(aes(y = 3), color = "#000", size = 1) +
  
  # SR line and ribbon
  geom_line(aes(y = est_SR), color = "#d62728", size = 1) +
 #eom_ribbon(aes(ymin = pci_lb_SR, ymax = pci_ub_SR), alpha = 0.2, fill = "#d62728") +
  
  labs(x = expression(delta),
       y = "Pooled Estimate",
       title = "Pooled Estimates with 95% CIs for CR and SR")+
  theme_minimal()

