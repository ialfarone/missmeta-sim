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
genimp.mnar = function(df, distribution = c("uniform", "normal", "tmvn"),
                  iter, minCR, maxCR, minSR, maxSR, meanCR, meanSR, sdCR, sdSR,
                  meantmv, sigmatmv, lower, upper, imprho, scaleSE) {
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
      cov = mv$vcov[1, 2],
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

true1 = mean(b3 + RTher[, 1])
true2 = mean(b3 + RTher[, 2])

results = rbind(
  sum.meth(resuni, true1, true2, "Uniform"),
  sum.meth(resnorm, true1, true2, "Normal(0,0)"),
  sum.meth(resnorm2, true1, true2, "Normal(3,3)"),
  sum.meth(resnorm3, true1, true2, "Normal(-3,-3)"),
  sum.meth(restmvn, true1, true2, "TMVN")
)

print(results)

resuni$method = "Uniform"
resnorm$method = "Normal(0,0)"
resnorm2$method = "Normal(3,3)"
resnorm3$method = "Normal(-3,-3)"
restmvn$method = "TMVN"

results.table = rbind(resuni, resnorm, resnorm2, resnorm3, restmvn)

# CR bias distribution
ggplot(results.table, aes(x = eff1 - true1, fill = method)) +
  geom_density(alpha = 0.4) +
  labs(title = "Bias Distribution for CR", x = "Bias", y = "Density") +
  theme_minimal()

# CR coverage 
library(tidyr)

resultslong = results %>%
  select(method, cover_CR, cover_SR) %>%
  pivot_longer(cols = c(cover_CR, cover_SR), names_to = "Outcome", values_to = "Coverage")

ggplot(resultslong, aes(x = method, y = Coverage, fill = Outcome)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Coverage Probability by Method", x = "Method", y = "Coverage") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme_minimal()

ggplot(results, aes(x = method)) +
  geom_linerange(aes(ymin = pci_lb_CR, ymax = pci_ub_CR), linewidth = 1.5) +
  geom_point(aes(y = est_CR), color = "black", size = 3) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Uncertainty Intervals for CR",
    x = "Method", y = "CR Interval"
  ) +
  theme_minimal()
