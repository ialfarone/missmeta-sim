# I simulate here a MNAR scenario for a mvmeta analysis of CR and SR.

# Cond 1: Missingness -> 5%, 10%, 15%, 20%
# Cond 2: Delta adjustment CR -> - 3, - 1.5, 0, 1.5, 3
# Cond 2bis: Delta adjustment SR -> - 3, - 1.5, 0, 1.5, 3
# Cond 3: Distribution -> uniform, normal(0), normal(3), normal(-3), mvt(3,3)

# I have then 4*5*5*5 = 500 situations

# I want to extract a list in which I have a dataset for each condition
# The goal is to look at the distribution of the CR and SR

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


genimp.mnar = function(df, distribution = c("uniform", "normal", "tmvn"),
                       iter, minCR, maxCR, minSR, maxSR, meanCR, meanSR, sdCR, sdSR,
                       meantmv, sigmatmv, lower, upper, imprho, delta, scaleSE) {
  distribution = match.arg(distribution)
  results = list()
  
  cor2cov = function(sd1, sd2, rho) {
    sd1 * sd2 * rho
  }
  
  for (d in delta) {
    for (i in 1:iter) {
      dfi = df 
      
      NmissCR = sum(is.na(dfi$EstCR))
      NmissSR = sum(is.na(dfi$EstSR))
      
      if (distribution == "uniform") {
        
        impCR = runif(NmissCR, min = minCR, max = maxCR)
        impSR = runif(NmissSR, min = minSR, max = maxSR)
        
      } else if (distribution == "normal") {
        
        impCR = rnorm(NmissCR, mean = meanCR + delta, sd = sdCR)
        impSR = rnorm(NmissSR, mean = meanSR, sd = sdSR)
        
      } else if (distribution == "tmvn") {
        
        imputed = rtmvnorm(
          n = sum(is.na(dfi$Cor.ws)),
          mean = meantmv + c(d, 0),
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
        delta = d,
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
