sim = function(seed = NULL, delta_seq = delta_seq,
                          target = target, meanCR, meanSR, sdCR, sdSR, 
                          iter = iter, true1 = 3, true2 = 3) {
  if (!is.null(seed)) set.seed(seed)
  
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
  
  # simulate ipd
  
  for (i in 1:S) {
    minA = runif(1, min = 18, max = 65)
    maxA = runif(1, min = minA + 5, max = 80)
    Age = runif(N, min = minA, max = maxA)
    pFem = 0.45
    Sex = rbinom(N, 1, pFem)
    Ther = rbinom(N, 1, 0.5)
    
    CR = rnorm(N, b0 + b1 * Age + b2 * Sex + (b3 + RTher[i, 1]) * Ther, sd = 5)
    SR = rnorm(N, b0 + b1 * Age + b2 * Sex + (b3 + RTher[i, 2]) * Ther, sd = 7)
    
    data[[i]] = data.frame(
      Study = i, Age = Age, Sex = factor(Sex), Ther = factor(Ther), CR = CR, SR = SR
    )
  }
  
  d = do.call(rbind, data)
  
  # pool estimates
  
  dat = lapply(1:S, function(s) {
    Sn = d[d$Study == s, ]
    fitsur = systemfit(list(CR = CR ~ Age + Sex + Ther, SR = SR ~ Age + Sex + Ther), 
                       "SUR", data = Sn)
    sum = summary(fitsur)
    data.frame(
      Study = s,
      EstCR = sum$coefficients[4, 1],
      EstSR = sum$coefficients[8, 1],
      SECR = sum$coefficients[4, 2],
      SESR = sum$coefficients[8, 2],
      Cor.ws = sum$residCor["CR", "SR"]
    )
  })
  
  dat = do.call(rbind, dat)
  
  # MNAR mechanism 
  
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
  meanCR = mean(dat$EstCR)
  meanSR = mean(dat$EstSR)
  
  for (i in conflict) {
    dist_cr = abs(dat$EstCR[i] - meanCR)
    dist_sr = abs(dat$EstSR[i] - meanSR)
    if (dist_cr > dist_sr) {
      M_sr[i] = 0  
    } else {
      M_cr[i] <- 0  
    }
  }
  
  dmnar = dat
  dmnar[M_cr == 1, c("EstCR", "SECR")] = NA
  dmnar[M_sr == 1, c("EstSR", "SESR")] = NA
  dmnar$Cor.ws[is.na(dmnar$EstCR) | is.na(dmnar$EstSR)] = NA
  
  # multiple imputation from donor distribution
   
  resnorm = genimp.mnar(
    df = dmnar,
    distribution = "normal",
    iter = iter,
    meanCR = meanCR,
    meanSR = meanSR,
    sdCR = sdCR,
    sdSR = sdSR,
    imprho = 0.7,
    delta = delta_seq,
    scaleSE = 1.5
  )
  
  # results by delta
  
  sumbydelta <- lapply(
    split(resnorm, resnorm$deltadj),
    function(df) sum.meth(df, true1 = true1, true2 = true2,
                          method_name = paste0("delta_", unique(df$deltadj)))
  )
  
  resbydelta = do.call(rbind, sumbydelta)
  resbydelta$replicate = seed
  return(resbydelta)
}

n_replicates <- 50

results_all <- lapply(1:n_replicates, function(seed) simulate_once(seed = seed))

results_combined <- do.call(rbind, results_all)

library(dplyr)
library(ggplot2)
library(tidyr)

results_combined <- results_combined %>%
  mutate(deltadj = as.numeric(gsub("delta_", "", method)))
ggplot(results_combined, aes(x = deltadj, y = est_CR, group = replicate, color = factor(replicate))) +
  geom_line(alpha = 0.6) +
  geom_point(alpha = 0.6) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.2, color = "black") +
  labs(x = expression(delta), y = "Estimated CR", color = "Replicate",
       title = "Estimated CR across delta (black = average)") +
  theme_minimal()

ggplot(results_combined, aes(x = deltadj, y = est_SR, group = replicate, color = factor(replicate))) +
  geom_line(alpha = 0.6) +
  geom_point(alpha = 0.6) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.2, color = "black") +
  labs(x = expression(delta), y = "Estimated SR", color = "Replicate",
       title = "Estimated SR across delta (black = average)") +
  theme_minimal()

bias_long <- results_combined %>%
  select(deltadj, replicate, bias_CR, bias_SR) %>%
  pivot_longer(cols = starts_with("bias"), names_to = "outcome", values_to = "bias") %>%
  mutate(outcome = recode(outcome, bias_CR = "CR", bias_SR = "SR"))

ggplot(bias_long, aes(x = deltadj, y = bias, color = outcome)) +
  geom_line(aes(group = interaction(replicate, outcome)), alpha = 0.5) +
  stat_summary(aes(group = outcome), fun = mean, geom = "line", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = expression(delta), y = "Bias", title = "Bias of Estimates across Delta and Replicates") +
  theme_minimal()

coverage_long <- results_combined %>%
  select(deltadj, replicate, cover_CR, cover_SR) %>%
  pivot_longer(cols = starts_with("cover"), names_to = "outcome", values_to = "coverage") %>%
  mutate(outcome = recode(outcome, cover_CR = "CR", cover_SR = "SR"))

ggplot(coverage_long, aes(x = deltadj, y = coverage, color = outcome)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.4) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen") +
  labs(x = expression(delta), y = "Coverage Probability", title = "Coverage Across Delta and Replicates") +
  theme_minimal()

