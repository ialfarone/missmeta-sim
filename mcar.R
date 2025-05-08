library(mvtnorm)
library(mixmeta)
library(systemfit)
library(dplyr)
library(tmvtnorm)

S = 1000
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

## Missing Completely At Random ####

mrate = 0.1
size = S * mrate
M_cr = sample(S, size = size, replace = T)
M_sr = sample(setdiff(1:S, M_cr), size = size, replace = FALSE)

dmcar = dat %>%
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
genimp = function(df,
                  distribution = c("uniform", "normal","tmvn"),
                  iter = NULL,
                  minCR = NULL,
                  maxCR = NULL,
                  minSR = NULL,
                  maxSR = NULL,
                  meanCR = NULL,
                  meanSR = NULL,
                  sdCR = NULL,
                  sdSR = NULL,
                  meantmv = NULL, 
                  sigmatmv = NULL, 
                  lower = NULL, 
                  upper = NULL,
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
    
    impSECR = rnorm(NmissCR, mean(dfi$SECR, na.rm = T), sd(dfi$SECR, na.rm = T)) 
    impSESR = rnorm(NmissSR, mean(dfi$SESR, na.rm = T), sd(dfi$SESR, na.rm = T))
    
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

# Occhio che https://cran.r-project.org/web/packages/metavcov/metavcov.pdf 
# incorpora gia' MI con multivariate meta analisi. Questo puo' dare un po'
# di solidita' e evidence.
# See also: Assessing the Sensitivity of Meta-analysis to Selection Bias: A Multiple Imputation
# Approach Author(s): James Carpenter, Gerta Rücker and Guido Schwarzer

#################################
#### Multivariate with MCAR #####
#################################

theta.m = cbind(dmcar$EstCR, dmcar$EstSR)
Sigma.m = cbind(dmcar$SECR^2,
              cor2cov(dmcar$SECR, dmcar$SESR, dmcar$Cor.ws),
              dmcar$SESR^2)

mv.m = mixmeta(theta.m, Sigma.m, method = "reml")
summary(mv.m)

# I wouldn't bother imputing

##################
#### Uniform #####
##################

NmissCR = sum(is.na(dmcar$EstCR))
NmissSR = sum(is.na(dmcar$EstSR))
impCR = runif(NmissCR, min = -max(d$CR), max = max(d$CR))
impSR = runif(NmissSR, min = -max(d$SR), max = max(d$SR))

plot(impCR)
plot(impSR)

resuni = genimp(
  df = dmcar,
  distribution = "uniform",
  iter = 10,
  minCR = -max(d$CR),
  maxCR = max(d$CR),
  minSR = -max(d$SR),
  maxSR = max(d$SR),
  imprho = 0.7
)

resuni

## Bias and Coverage from Rubin's rules

### Rubin (IA: write formulae somewhere)
m = nrow(resuni)

peff1 = mean(resuni$eff1) # pooled effects
peff2 = mean(resuni$eff2)

pse1 = mean(resuni$se1^2)
pse2 = mean(resuni$se2^2)

btwvar1 = var(resuni$eff1)
btwvar2 = var(resuni$eff2)

totvar1 = pse1 + (1 + 1/m) * btwvar1
totvar2 = pse2 + (1 + 1/m) * btwvar2

pse1 = sqrt(totvar1) # se pooled
pse2 = sqrt(totvar2)

pci1 = peff1 + c(-1, 1) * qnorm(0.975) * pse1
pci2 = peff2 + c(-1, 1) * qnorm(0.975) * pse2

#### bias and coverage

true1 = mean(b3 + RTher[, 1])
true2 = mean(b3 + RTher[, 2])

bias1 = mean(resuni$eff1) - true1
bias2 = mean(resuni$eff2) - true2

cover1 = mean(resuni$ci.lb1 <= true1 & resuni$ci.ub1 >= true1)
cover2 = mean(resuni$ci.lb2 <= true2 & resuni$ci.ub2 >= true2)

bias1
bias2

cover1
cover2

pci1
pci2

hist(resuni$eff1 - true1, breaks = 50, main = "Bias Distribution (CR)", xlab = "Bias")
abline(v = bias1, col = "red", lwd = 2)

################
#### Normal ####
################

meanCR = 0
meanSR = 0
sdCR = 10
sdSR = 12

impCR = rnorm(NmissCR, mean = meanCR, sd = sdCR)
impSR = rnorm(NmissSR, mean = meanSR, sd = sdSR)
plot(impCR)

resnorm = genimp(
  df = dmcar,
  distribution = "normal",
  iter = 20,
  meanCR = 0, meanSR = 0,
  sdCR = 10, sdSR = 12,
  imprho = 0.7
)
resnorm

## Bias and Coverage from Rubin's rules

### Rubin (IA: write formulae somewhere)
m = nrow(resnorm)

peff1 = mean(resnorm$eff1) # pooled effects
peff2 = mean(resnorm$eff2)

pse1 = mean(resnorm$se1^2)
pse2 = mean(resnorm$se2^2)

btwvar1 = var(resnorm$eff1)
btwvar2 = var(resnorm$eff2)

totvar1 = pse1 + (1 + 1/m) * btwvar1
totvar2 = pse2 + (1 + 1/m) * btwvar2

pse1 = sqrt(totvar1) # se pooled
pse2 = sqrt(totvar2)

pci1 = peff1 + c(-1, 1) * qnorm(0.975) * pse1
pci2 = peff2 + c(-1, 1) * qnorm(0.975) * pse2

#### bias and coverage

true1 = mean(b3 + RTher[, 1])
true2 = mean(b3 + RTher[, 2])

bias1 = mean(resnorm$eff1) - true1
bias2 = mean(resnorm$eff2) - true2

cover1 = mean(resnorm$ci.lb1 <= true1 & resnorm$ci.ub1 >= true1)
cover2 = mean(resnorm$ci.lb2 <= true2 & resnorm$ci.ub2 >= true2)

bias1
bias2

cover1
cover2

pci1
pci2

hist(resnorm$eff1 - true1, breaks = 50, main = "Bias Distribution (CR)", xlab = "Bias")
abline(v = bias1, col = "red", lwd = 2)


resnorm2 = genimp(
  df = dmcar,
  distribution = "normal",
  iter = 20,
  meanCR = 3, meanSR = 3,
  sdCR = 10, sdSR = 12,
  #  impSECR = 100,
  #  impSESR = 100,
  imprho = 0.7
)
resnorm2

resnorm3 = genimp(
  df = dmcar,
  distribution = "normal",
  iter = 20,
  meanCR = -3, meanSR = -3,
  sdCR = 10, sdSR = 12,
  imprho = 0.7
)
resnorm3

#############################
##### Truncate MVNormal #####
#############################

lower = c(-Inf, -Inf) # for now simply multivariate
upper = c(Inf, Inf)

# for now from the data, but this does not hold for MNAR or MAR 

meantmv = colMeans(cbind(dmcar$EstCR, dmcar$EstSR), na.rm = TRUE)
sigmatmv = cov(cbind(dmcar$EstCR, dmcar$EstSR), use = "complete.obs")

restmvn = genimp(
  df = dmcar,
  distribution = "tmvn",
  iter = 10,
  lower = lower,
  upper = upper, 
  meantmv = meantmv,
  sigmatmv = sigmatmv, 
  imprho = 0.7
)

restmvn
