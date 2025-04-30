library(mvtnorm)
library(mixmeta)
library(systemfit)
library(dplyr)

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

## Missing At Random ####

sub = d %>%
  group_by(Study) %>%
  summarise(meanA = mean(Age))

beta0 = qlogis(0.10)
beta1 = 0.15  # higher age -> higher prob of CR being observed

plot(sub$meanA, plogis(beta0 + beta1 * (sub$meanA - 18)))

# sub$likObs = 1 / (1 + exp(-(beta0 + beta1 * sub$meanA)))

# centering on the minimum thus 18 is 0
sub$likObs <- plogis(beta0 + beta1 * (sub$meanA - 18))

M0 = rbinom(S, size = 1, prob = sub$likObs)
dmar = dat

for (i in 1:S) {
  if (M0[i] == 0) {
    dmar$EstCR[dat$Study == i] = NA
    dmar$SECR[dat$Study == i] = NA
    dmar$Cor.ws[dat$Study == i] = NA
    
  }
}

head(dmar)

ggplot(dmar, aes(x = sub$meanA, fill = is.na(EstCR))) +
  geom_histogram(binwidth = 2, position = "stack") +
  labs(title = "Missingness in CR as a function of Age", x = "Age", y = "Count") +
  scale_fill_manual(values = c("#5B5F97", "#EE6C4D"), name = "CR Missing") +
  theme_minimal()

# Random sample generator for continuous missing data ####
genimp = function(df,
                  distribution = c("uniform", "normal"),
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
                  imprho = NULL) {
  distribution = match.arg(distribution)
  results = vector(mode = "list", length = iter)
  
  cor2cov = function(sd1, sd2, rho) {
    sd1 * sd2 * rho
  }
  
  for (i in 1:iter) {
    dfi = df
    
    NmissCR = sum(is.na(dfi$EstCR))
    NmissSR = sum(is.na(dfi$EstSR))
    
    # Draw imputed values based on selected distribution
    if (distribution == "uniform") {
      if (is.null(minCR) || is.null(maxCR) ||
          is.null(minSR) || is.null(maxSR)) {
        stop(
          "For uniform distribution, 'minCR', 'maxCR', 'minSR', and 'maxSR'
            must all be provided."
        )
      }
      
      impCR = runif(NmissCR, min = minCR, max = maxCR)
      impSR = runif(NmissSR, min = minSR, max = maxSR)
      
    } else if (distribution == "normal") {
      if (is.null(sdSR) || is.null(sdCR) ||
          is.null(meanSR) || is.null(meanCR)) {
        stop("For normal distribution, 'sdCR' and 'sdSR' must both be provided.")
      }
      
      impCR <- rnorm(NmissCR, mean = meanCR, sd = sdCR)
      impSR <- rnorm(NmissSR, mean = meanSR, sd = sdSR)
    }
    
    dfi$EstCR[is.na(dfi$EstCR)] = impCR
    dfi$EstSR[is.na(dfi$EstSR)] = impSR
    #    dmcar$SECR[is.na(dmcar$SECR)] = impSECR
    #    dmcar$SESR[is.na(dmcar$SESR)] = impSESR
    dfi$Cor.ws[is.na(dfi$Cor.ws)] = imprho
    
    impSECR = sample(dfi$SECR[!is.na(dfi$SECR)], NmissCR, replace = TRUE)
    impSESR = sample(dfi$SESR[!is.na(dfi$SESR)], NmissSR, replace = TRUE)
    
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

resuni = genimp(
  df = dmar,
  distribution = "uniform",
  iter = 10,
  minCR = -max(d$CR),
  maxCR = max(d$CR),
  minSR = -max(d$SR),
  maxSR = max(d$SR),
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7
)
resuni

resnorm = genimp(
  df = dmar,
  distribution = "normal",
  iter = 20,
  meanCR = 0,
  meanSR = 0,
  sdCR = 10,
  sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7
)
resnorm

resnorm2 = genimp(
  df = dmar,
  distribution = "normal",
  iter = 20,
  meanCR = 3,
  meanSR = 3,
  sdCR = 10,
  sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7
)
resnorm2

resnorm3 = genimp(
  df = dmcar,
  distribution = "normal",
  iter = 20,
  meanCR = -3,
  meanSR = -3,
  sdCR = 10,
  sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7
)
resnorm3
