##################
rm(list=ls())
##################

library(mvtnorm)
library(mixmeta)
library(systemfit)
library(dplyr)
library(ggplot2)
library(tmvtnorm)
library(missmeta)


# Simulate dataset ####

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
  maxA = runif(1, min = minA + 5, max = 80)
  Age = runif(N, min = minA, max = maxA)
  
  pFem = 0.45
  Sex = rbinom(N, size = 1, prob = pFem)
  Ther = rbinom(N, size = 1, prob = 0.5)
  
  CR = rnorm(N,
             mean = b0 + b1 * Age + b2 * Sex + (b3 + RTher[i, 1]) * Ther,
             sd = 5)
  SR = rnorm(N,
             mean = b0 + b1 * Age + b2 * Sex + (b3 + RTher[i, 2]) * Ther,
             sd = 7)
  
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

# Calculate effect size and standard errors ####

for (s in 1:S) {
  Sn = d[d$Study == s, ]
  
  m1 = CR ~ Age + Sex + Ther
  m2 = SR ~ Age + Sex + Ther
  
  fitsur = systemfit(list(CR = m1, SR = m2), "SUR", data = Sn)
  sum = summary(fitsur)
  
  dat[[s]] = data.frame(
    Study = s,
    N = N,
    EstCR = sum$coefficients[4, 1],
    EstSR = sum$coefficients[8, 1],
    SECR = sum$coefficients[4, 2],
    SESR = sum$coefficients[8, 2],
    Cor.ws = sum$residCor["CR", "SR"]
  )
}

dat = do.call(rbind, dat)
head(dat)

# Multivariate meta-analysis ####

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
betaCR = 2
betaSR = -1.5
invlogit <- function(x) plogis(x)

beta0_cr = uniroot(function(b0) mean(invlogit(b0 + betaCR * dat$EstCR)) - (1 - target / 2), c(-20, 20))$root
prob_cr = invlogit(beta0_cr + betaCR * dat$EstCR)
beta0_sr = uniroot(function(b0) mean(invlogit(b0 + betaSR * dat$EstSR)) - (1 - target / 2), c(-20, 20))$root
prob_sr = invlogit(beta0_sr + betaSR * dat$EstSR)

M_cr = rbinom(nrow(dat), 1, 1 - prob_cr)
M_sr = rbinom(nrow(dat), 1, 1 - prob_sr)

conflict = which(M_cr == 1 & M_sr == 1)
meanCR_obs = mean(dat$EstCR)
meanSR_obs = mean(dat$EstSR)

for (i in conflict) {
  dist_cr = abs(dat$EstCR[i] - meanCR_obs)
  dist_sr = abs(dat$EstSR[i] - meanSR_obs)
  if (dist_cr > dist_sr) {
    M_sr[i] = 0
  } else {
    M_cr[i] = 0
  }
}

dmnar = dat

dmnar[M_cr == 1, c("EstCR", "SECR")] = NA
dmnar[M_sr == 1, c("EstSR", "SESR")] = NA
dmnar$Cor.ws[is.na(dmnar$EstCR) | is.na(dmnar$EstSR)] = NA

head(dmnar)
sum(is.na(dmnar$Cor.ws))


###############
# genimp mod
###############

genimp = function(df, iter = 100,
                   imp1, imp2,
                   eff1 = "Eff1", eff2 = "Eff2",
                   se1 = "SE1", se2 = "SE2",
                   cor = "Cor.ws", N = "N",
                   imprho = 0.7) {
  results = vector("list", iter)
  
  for (i in seq_len(iter)) {
    dfi = df
    Nmiss1 = sum(is.na(dfi[[eff1]]))
    Nmiss2 = sum(is.na(dfi[[eff2]]))
    
    dfi[[eff1]][is.na(dfi[[eff1]])] = imp1(Nmiss1)
    dfi[[eff2]][is.na(dfi[[eff2]])] = imp2(Nmiss2)
    dfi[[cor]][is.na(dfi[[cor]])] = imprho
    
    comp = !is.na(dfi[[se1]]) & !is.na(dfi[[se2]])
    log_sig = cbind(
      log(dfi[[se1]][comp] * sqrt(dfi[[N]][comp])),
      log(dfi[[se2]][comp] * sqrt(dfi[[N]][comp]))
    )
    
    
    mu_hat    = colMeans(log_sig)
    Sigma_hat = var(log_sig)     
    
    miss = is.na(dfi[[se1]]) | is.na(dfi[[se2]])
    log_draw  = mvtnorm::rmvnorm(sum(miss), mean = mu_hat, sigma = Sigma_hat)
    sigma_imp = exp(log_draw)    
    
    idx_se1 = which(miss & is.na(dfi[[se1]]))
    idx_se2 = which(miss & is.na(dfi[[se2]]))
    
    dfi[[se1]][idx_se1] = sigma_imp[seq_along(idx_se1), 1] / sqrt(dfi[[N]][idx_se1])
    dfi[[se2]][idx_se2] = sigma_imp[seq_along(idx_se2), 2] / sqrt(dfi[[N]][idx_se2])
    
    results[[i]] <- dfi
  }
  
  results
}


imp1 <- function(n) rnorm(n, mean = 3, sd = 2)
imp2 <- function(n) runif(n, min = 4, max = 7)

sigma = matrix(c(6^2,
                 CorCov(6,6,0.7), 
                 CorCov(6,6,0.7), 
                 6^2), nrow = 2)

lower = c(-100, -100)
upper = c(100, 100)

imp1 = function(n) rtmvnorm(n, mean = c(1, 6), sigma = sigma, lower = lower, upper = upper)[,1]
imp2 = function(n) rtmvnorm(n, mean = c(1, 6), sigma = sigma, lower = lower, upper = upper)[,2]

out <- genimp(df = dmnar,
              iter = 2,
              imp1 = imp1,
              imp2 = imp2,
              eff1 = "EstCR",
              eff2 = "EstSR",
              se1 = "SECR",
              se2 = "SESR",
              cor = "Cor.ws",
              N = "N",
              imprho = 0.3)
out
