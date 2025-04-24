library(mvtnorm)
library(mixmeta)

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


dat = list()

for (i in 1:S) {
  
  minA = runif(1, min = 18, max = 65)
  maxA = runif(1, min = minA + 5, max = 90)
  Age = runif(N, min = minA, max = maxA)
  
  pFem = 0.45
  Sex = rbinom(N, size = 1, prob = pFem)  
  Ther = rbinom(N, size = 1, prob = 0.5)  
  
  CR = rnorm(N, mean = b0 + b1*Age + b2*Sex + (b3 + RTher[i, 1])*Ther, sd = 5)
  SR = rnorm(N, mean = b0 + b1*Age + b2*Sex + (b3 + RTher[i, 2])*Ther, sd = 9)
  
  dat[[i]] = data.frame(
    Study = i,
    Age = Age,
    Sex = factor(Sex, levels = 0:1, labels = c("M", "F")),
    Ther = factor(Ther, levels = 0:1, labels = c("New", "Std")),
    CR = CR,
    SR = SR
  )
}

d = do.call(rbind, dat)

library(systemfit)

data <- data.frame(
  CR_Therapy = numeric(S),
  SR_Therapy = numeric(S),
  theta_1 = numeric(S),
  theta_2 = numeric(S),
  theta_12 = numeric(S),
  cor = numeric(S)
)

for (s in 1:S) {
  Sn <- d[d$Study == s, ]
  
  m1 <- CR ~ Age + Sex + Ther 
  m2 <- SR ~ Age + Sex + Ther
  
  fitsur <- systemfit(list(CR = m1, SR = m2), "SUR", data = Sn)
  sum <- summary(fitsur)
  
  data$CR_Therapy[s] <- sum$coefficients[4, 1]
  data$SR_Therapy[s] <- sum$coefficients[8, 1]
  
  data$theta_1[s] <- sum$coefficients[4, 2]
  data$theta_2[s] <- sum$coefficients[8, 2]
  
  data$cor[s] <- sum$residCor["CR", "SR"]
}

dat = as.data.frame(cbind(1:S, N, data$CR_Therapy, data$SR_Therapy, 
                          data$theta_1, data$theta_2, data$cor))

colnames(dat) <- c("Study", "N", "EstCR", "EstSR", "SECR", "SESR", "Cor.ws")
head(dat)

## Multivariate meta-analysis ####

theta <- cbind(dat$EstCR, dat$EstSR)
cor2cov = function(sd1, sd2, rho) {sd1 * sd2 * rho}
Sigma = cbind(dat$SECR^2, cor2cov(dat$SECR, dat$SESR, dat$Cor.ws), dat$SESR^2)

mv.c <- mixmeta(theta, Sigma, method="reml")
summary(mv.c)

# Generate missing ####

## Missing Completely At Random ####

mrate = 0.2
size = S*mrate
M_cr = sample(S, size=size, replace=T)
M_sr = ifelse(S != M_cr, sample(S, size=size, replace=T))

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


## Missing At Random ####

sub <- d %>%
  group_by(Study) %>%
  summarise(meanA = mean(Age))

beta0 = 10     
beta1 = -0.15  # higher age -> lower prob of CR being observed

sub$likObs = 1 / (1 + exp(-(beta0 + beta1 * sub$meanA)))

M0 <- rbinom(100, size = 1, prob = sub$likObs)
dmar = dat

for (i in 1:S) {
  if (M0[i] == 0) {
    dmar$EstCR[dat$Study == i] <- NA
  }
}
head(dmar)

library(ggplot2)

ggplot(dmar, aes(x = sub$meanA, fill = is.na(EstCR))) +
  geom_histogram(binwidth = 2, position = "stack") +
  labs(title = "Missingness in CR as a function of Age", x = "Age", y = "Count") +
  scale_fill_manual(values = c("#5B5F97", "#EE6C4D"), name = "CR Missing") +
  theme_minimal()


## Missing Not At Random (focused) ####

beta0 = 0     
beta1 = 0.30  # lower EstCR -> lower prob of CR being observed

dmnar = dat

likObs = 1 / (1 + exp(-(beta0 + beta1 * dmnar$EstCR)))

M0 <- rbinom(100, size = 1, prob = likObs)

for (i in 1:S) {
  if (M0[i] == 0) {
    dmnar$EstCR[dat$Study == i] <- NA
  }
}
head(dmnar)
sum(is.na(dmnar$EstCR))

library(ggplot2)

ggplot(dmnar, aes(x = dat$EstCR, fill = is.na(EstCR))) +
  geom_histogram(binwidth = 2, position = "stack") +
  labs(title = "Missingness in CR as a function of CR", x = "EstCR", y = "Count") +
  scale_fill_manual(values = c("#5B5F97", "#EE6C4D"), name = "CR Missing") +
  theme_minimal()

# Random sample generator for continuous missing data ####

## Uniform distribution ####

results <- list()
iter <- 10

dmcar_orig <- dmcar 

for (i in 1:iter) {

  dmcar = dmcar_orig
  
  unif_cr = runif(sum(is.na(dmcar$EstCR)), min = -max(d$CR), max = max(d$CR))
  unif_sr = runif(sum(is.na(dmcar$EstSR)), min = -max(d$SR), max = max(d$SR))
  
  dmcar$EstCR[is.na(dmcar$EstCR)] = unif_cr
  dmcar$EstSR[is.na(dmcar$EstSR)] = unif_sr
  
  theta = cbind(dmcar$EstCR, dmcar$EstSR)
  
  Sigma = cbind(dmcar$SECR^2, cor2cov(dmcar$SECR, dmcar$SESR, dmcar$Cor.ws), 
                 dmcar$SESR^2)

  
  mv.c = mixmeta(theta, Sigma, method = "reml")
  
  results[[i]] = mv.c$coefficients
}

do.call(rbind, results)

