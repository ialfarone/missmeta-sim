library(mvtnorm)
library(mixmeta)
library(systemfit)
library(dplyr)
library(ggplot2)
library(tmvtnorm)

S = 100
N = 1000

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
  Sn <- d[d$Study == s, ]
  
  m1 <- CR ~ Age + Sex + Ther
  m2 <- SR ~ Age + Sex + Ther
  
  fitsur <- systemfit(list(CR = m1, SR = m2), "SUR", data = Sn)
  sum <- summary(fitsur)
  
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

target = 0.2  
beta1 = 2

beta0 = uniroot(function(b0) mean(plogis(b0 + beta1 * dat$EstCR)) -
                   (1 - target),
                 c(-20, 20))$root

prob_obs = plogis(beta0 + beta1 * dat$EstCR)

dmnar    = dat
M0       = rbinom(nrow(dat), 1, prob_obs)
dmnar[M0 == 0, c("EstCR", "SECR", "Cor.ws")] = NA

head(dmnar)
sum(is.na(dmnar$EstCR))

ggplot(dmnar, aes(x = dat$EstCR, fill = is.na(EstCR))) +
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


# Random sample generator for continuous missing data for MNAR ####
genimp.mnar = function(df,
                  distribution = c("uniform", "normal", "tmvn"),
                  iter = NULL,
                  minCR = NULL,
                  maxCR = NULL,
                  minSR = NULL,
                  maxSR = NULL,
                  meanCR = NULL,
                  meanSR = NULL,
                  sdCR = NULL,
                  sdSR = NULL,
                  #                  impSECR = NULL,
                  #                  impSESR = NULL,
                  meantmv = NULL, 
                  sigmatmv = NULL, 
                  lower = lower, 
                  upper = upper, 
                  imprho = NULL,
                  scaleSE = NULL) {
  distribution = match.arg(distribution)
  results = vector(mode = "list", length = iter)
  
  cor2cov = function(sd1, sd2, rho) {
    sd1 * sd2 * rho
  }
  
  for (i in 1:iter) {
    dfi = df 
    
    NmissCR = sum(is.na(dfi$EstCR))
    NmissSR = sum(is.na(dfi$EstSR))
    
    if (distribution == "uniform") {
      
      impCR = runif(NmissCR, min = minCR, max = maxCR)
      impSR = runif(NmissSR, min = minSR, max = maxSR)
      
    } else if (distribution == "normal") {
      
      impCR = rnorm(NmissCR, mean = meanCR, sd = sdCR)
      impSR = rnorm(NmissSR, mean = meanSR, sd = sdSR)
      
    } else if (distribution == "tmvn") {
      
      imputed = rtmvnorm(
        n = sum(is.na(dfi$Cor.ws)),
        mean = meantmv,
        sigma = sigmatmv,
        lower = lower,
        upper = upper
      )
      
      impCR = imputed[, 1]
      impSR = imputed[, 2]
    }
    
    
    dfi$EstCR[is.na(dfi$EstCR)] = sample(impCR, NmissCR, replace = F)
    dfi$EstSR[is.na(dfi$EstSR)] = sample(impSR, NmissSR, replace = F)
    #    dmcar$SECR[is.na(dmcar$SECR)] = impSECR
    #    dmcar$SESR[is.na(dmcar$SESR)] = impSESR
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
    
    results[[i]] = data.frame(
      iter = i,
      eff1 = mv$coefficients[1],
      eff2 = mv$coefficients[2],
      se1 = sqrt(mv$vcov[1, 1]),
      se2 = sqrt(mv$vcov[2, 2]),
      ci.lb1 = ci[1, 1],
      ci.ub1 = ci[1, 2],
      ci.lb2 = ci[2, 1],
      ci.ub2 = ci[2, 2]
    )
  }
  
  do.call(rbind, results)
}

resuni = genimp.mnar(
  df = dmnar,
  distribution = "uniform",
  iter = 10,
  minCR = -max(d$CR),
  maxCR = max(d$CR),
  minSR = -max(d$SR),
  maxSR = max(d$SR),
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7,
  scaleSE = 1.5
)
resuni

resnorm = genimp.mnar(
  df = dmnar,
  distribution = "normal",
  iter = 10,
  meanCR = 0,
  meanSR = 0,
  sdCR = 10,
  sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7,
  scaleSE = 1.5
)
resnorm

resnorm2 = genimp.mnar(
  df = dmnar,
  distribution = "normal",
  iter = 10,
  meanCR = 3,
  meanSR = 3,
  sdCR = 10,
  sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7,
  scaleSE = 1.5
)
resnorm2

resnorm3 = genimp.mnar(
  df = dmnar,
  distribution = "normal",
  iter = 10,
  meanCR = -3,
  meanSR = -3,
  sdCR = 10,
  sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7,
  scaleSE = 1.5
)
resnorm3

#############################
##### Truncate MVNormal #####
#############################
lower = c(-Inf, -Inf) # for now simply multivariate
upper = c(Inf, Inf)

# for now from the data, but this does not hold for MNAR  

meantmv = colMeans(cbind(dmnar$EstCR, dmnar$EstSR), na.rm = TRUE)
sigmatmv = cov(cbind(dmnar$EstCR, dmnar$EstSR), use = "complete.obs")

restmvn = genimp.mnar(
  df = dmnar,
  distribution = "tmvn",
  iter = 10,
  lower = lower,
  upper = upper, 
  meantmv = meantmv,
  sigmatmv = sigmatmv, 
  imprho = 0.7,
  scaleSE = 1.5
)

restmvn


################################################################################
# Calculate bias and coverage and compare

evaluate_method <- function(res, true1, true2, method_name) {
  m <- nrow(res)
  
  peff1 <- mean(res$eff1)
  peff2 <- mean(res$eff2)
  
  pse1 <- mean(res$se1^2)
  pse2 <- mean(res$se2^2)
  
  btwvar1 <- var(res$eff1)
  btwvar2 <- var(res$eff2)
  
  totvar1 <- pse1 + (1 + 1/m) * btwvar1
  totvar2 <- pse2 + (1 + 1/m) * btwvar2
  
  pse1 <- sqrt(totvar1)
  pse2 <- sqrt(totvar2)
  
  pci1 <- peff1 + c(-1, 1) * qnorm(0.975) * pse1
  pci2 <- peff2 + c(-1, 1) * qnorm(0.975) * pse2
  
  bias1 <- peff1 - true1
  bias2 <- peff2 - true2
  
  cover1 <- mean(res$ci.lb1 <= true1 & res$ci.ub1 >= true1)
  cover2 <- mean(res$ci.lb2 <= true2 & res$ci.ub2 >= true2)
  
  return(data.frame(
    method = method_name,
    bias_CR = bias1,
    bias_SR = bias2,
    cover_CR = cover1,
    cover_SR = cover2,
    pci_lb_CR = pci1[1],
    pci_ub_CR = pci1[2],
    pci_lb_SR = pci2[1],
    pci_ub_SR = pci2[2]
  ))
}

true1 <- mean(b3 + RTher[, 1])
true2 <- mean(b3 + RTher[, 2])

results_all <- rbind(
  evaluate_method(resuni, true1, true2, "Uniform"),
  evaluate_method(resnorm, true1, true2, "Normal(0,0)"),
  evaluate_method(resnorm2, true1, true2, "Normal(3,3)"),
  evaluate_method(resnorm3, true1, true2, "Normal(-3,-3)"),
  evaluate_method(restmvn, true1, true2, "TMVN")
)

print(results_all)

resuni$method <- "Uniform"
resnorm$method <- "Normal(0,0)"
resnorm2$method <- "Normal(3,3)"
resnorm3$method <- "Normal(-3,-3)"
restmvn$method <- "TMVN"

res_all <- rbind(resuni, resnorm, resnorm2, resnorm3, restmvn)

# CR bias distribution
ggplot(res_all, aes(x = eff1 - true1, fill = method)) +
  geom_density(alpha = 0.4) +
  labs(title = "Bias Distribution for CR", x = "Bias", y = "Density") +
  theme_minimal()

# CR coverage 

library(tidyr)

results_all_long <- results_all %>%
  select(method, cover_CR, cover_SR) %>%
  pivot_longer(cols = c(cover_CR, cover_SR), names_to = "Outcome", values_to = "Coverage")

ggplot(results_all_long, aes(x = method, y = Coverage, fill = Outcome)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Coverage Probability by Method", x = "Method", y = "Coverage") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme_minimal()


