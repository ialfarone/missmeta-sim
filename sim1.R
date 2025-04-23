library(mvtnorm)
library(mixmeta)

S = 500
N = 1000  

Mu = c(0, 0)
Tau = c(1, 2)
rho = 0.7
Sigma = diag(Tau) %*% matrix(c(1, rho, rho, 1), nrow = 2) %*% diag(Tau)
RTher = rmvnorm(S, Mu, Sigma)

b0 = 20
b1 = 0.5
b2 = 1.5
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
  SR = rnorm(N, mean = b0 + b1*Age + b2*Sex + (b3 + RTher[i, 2])*Ther, sd = 5)
  
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

dat = as.data.frame(cbind(N, data$CR_Therapy, data$SR_Therapy, 
                          data$theta_1, data$theta_2, data$cor))

colnames(dat) <- c("N", "EstCR", "EstSR", "SECR", "SESR", "Cor.ws")
head(dat)

theta <- cbind(dat$EstCR, dat$EstSR)
cor2cov = function(sd1, sd2, rho) {sd1 * sd2 * rho}
Sigma = cbind(dat$SECR^2, cor2cov(dat$SECR, dat$SESR, dat$Cor.ws), dat$SESR^2)

mv.c <- mixmeta(theta, Sigma, method="reml")
summary(mv.c)
