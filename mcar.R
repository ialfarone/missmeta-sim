library(mvtnorm)
library(mixmeta)
library(systemfit)

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

dat = vector(mode = "list", length = S)

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

data = vector(mode = "list", length = S)

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

theta <- cbind(dat$EstCR, dat$EstSR)
cor2cov = function(sd1, sd2, rho) {
  sd1 * sd2 * rho
}
Sigma = cbind(dat$SECR^2,
              cor2cov(dat$SECR, dat$SESR, dat$Cor.ws),
              dat$SESR^2)

mv.c <- mixmeta(theta, Sigma, method = "reml")
summary(mv.c)

# Generate missing #### 

## Missing Completely At Random ####

mrate = 0.1
size = S * mrate
M_cr = sample(S, size = size, replace = T)
M_sr = sample(setdiff(1:S, M_cr), size = size, replace = FALSE)

library(dplyr)
dmcar <- dat %>%
  mutate(
    EstCR = ifelse(Study %in% M_cr, NA, EstCR),
    SECR = ifelse(Study %in% M_cr, NA, SECR),
    EstSR = ifelse(Study %in% M_sr, NA, EstSR),
    SESR = ifelse(Study %in% M_sr, NA, SESR),
    Cor.ws = ifelse(Study %in% M_cr | Study %in% M_sr, NA, Cor.ws)
  )

dmcar

sum(is.na(dmcar$Cor.ws)) / S

library(ggplot2)

ggplot(dmcar, aes(x = dat$EstCR, fill = is.na(EstCR))) +
  geom_histogram(binwidth = 2, position = "stack") +
  labs(title = "Missingness in CR", x = "EstCR", y = "Count") +
  scale_fill_manual(values = c("#5B5F97", "#EE6C4D"), name = "CR Missing") +
  theme_minimal()

# Random sample generator for continuous missing data ####

## Uniform distribution ####

iter = 10
results = vector(mode = "list", length = iter)

dmcar_orig = dmcar

for (i in 1:iter) {
  dmcar = dmcar_orig
  
  unif_cr = runif(sum(is.na(dmcar$EstCR)), min = -max(d$CR), max = max(d$CR))
  unif_sr = runif(sum(is.na(dmcar$EstSR)), min = -max(d$SR), max = max(d$SR))
  
  dmcar$EstCR[is.na(dmcar$EstCR)] = unif_cr
  dmcar$EstSR[is.na(dmcar$EstSR)] = unif_sr
  
  # check this
  dmcar$SECR[is.na(dmcar$SECR)] = 10^2
  dmcar$SESR[is.na(dmcar$SESR)] = 10^2
  dmcar$Cor.ws[is.na(dmcar$Cor.ws)] = 0.7
  
  theta = cbind(dmcar$EstCR, dmcar$EstSR)
  
  Sigma = cbind(dmcar$SECR^2,
                cor2cov(dmcar$SECR, dmcar$SESR, dmcar$Cor.ws),
                dmcar$SESR^2)
  
  mv = mixmeta(theta, Sigma, method = "reml")
  ci = confint(mv)
  
  results[[i]] = data.frame(
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

res = do.call(rbind, results)
res
### Option 2: normal distribution with mean 0 and sd 10 and sd 12

results = vector(mode = "list", length = iter)

dmcar_orig = dmcar

for (i in 1:iter) {
  dmcar = dmcar_orig
  
  unif_cr = rnorm(sum(is.na(dmcar$EstCR)), mean = 0, sd = 10)
  unif_sr = rnorm(sum(is.na(dmcar$EstSR)), mean = 0, sd = 12)
  
  dmcar$EstCR[is.na(dmcar$EstCR)] = unif_cr
  dmcar$EstSR[is.na(dmcar$EstSR)] = unif_sr
  
  # check this
  dmcar$SECR[is.na(dmcar$SECR)] = 10^2
  dmcar$SESR[is.na(dmcar$SESR)] = 10^2
  dmcar$Cor.ws[is.na(dmcar$Cor.ws)] = 0.7
  
  theta = cbind(dmcar$EstCR, dmcar$EstSR)
  
  Sigma = cbind(dmcar$SECR^2,
                cor2cov(dmcar$SECR, dmcar$SESR, dmcar$Cor.ws),
                dmcar$SESR^2)
  
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

res2 = do.call(rbind, results)

hist(res2[, 'y1'])
hist(res2[, 'y2'])
