################
rm(list = ls())
################

library(mixmeta)
library(mvtnorm)
library(systemfit)  
library(tmvtnorm)   
library(ggplot2)
library(dplyr)

##############################
# genimp.mnar: Imputation ####
##############################

cor2cov = function(sd1, sd2, rho) sd1 * sd2 * rho

genimp.mnar = function(df, iter = 100, distribution = c("uniform", "normal", "tmvn"),
                       minCR, maxCR, minSR, maxSR, meanCR, meanSR, sdCR, sdSR,
                       meantmv, sigmatmv, lower, upper, imprho) {
  
  distribution = match.arg(distribution)
  results = vector("list", iter)
  
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
      
      dfi$EstCR[is.na(dfi$EstCR)] = impCR
      dfi$EstSR[is.na(dfi$EstSR)] = impSR
      dfi$Cor.ws[is.na(dfi$Cor.ws)] = imprho
      
      # Irene: See here from e-mail
      # new imputation of std.err from lognormal biv distr (see syst rev)
      
      log_sig = with(dfi[!is.na(dfi$SECR) & !is.na(dfi$SESR), ],
                      cbind(log(SECR * sqrt(N)),
                            log(SESR * sqrt(N))))
      mu_hat    = colMeans(log_sig)
      Sigma_hat = var(log_sig)     
      
      miss = with(dfi, is.na(SECR) | is.na(SESR))
      log_draw  = rmvnorm(sum(miss), mu_hat, Sigma_hat)
      sigma_imp = exp(log_draw)    
      
      idx_CR = miss & is.na(dfi$SECR)
      idx_SR = miss & is.na(dfi$SESR)
      
      dfi$SECR[idx_CR] = sigma_imp[ idx_CR[miss], 1 ] / sqrt(dfi$N[idx_CR])
      dfi$SESR[idx_SR] = sigma_imp[ idx_SR[miss], 2 ] / sqrt(dfi$N[idx_SR])
      
      # Irene: from here same as before
      
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
        ci.lb1 = ci[1, 1], ci.ub1 = ci[1, 2],
        ci.lb2 = ci[2, 1], ci.ub2 = ci[2, 2]
      )
    }
  do.call(rbind, results)
}

##############################
# sum.meth: Rubin's rules ####
##############################

sum.meth = function(res, true1, true2) {
  m = nrow(res)
  Q_mat = cbind(res$eff1, res$eff2)
  Q_bar = colMeans(Q_mat)
  B = cov(Q_mat)
  
  U_list = lapply(1:m, function(i) {
    matrix(c(res$se1[i]^2, res$cov[i], res$cov[i], res$se2[i]^2), nrow = 2)
  })
  U_bar = Reduce("+", U_list) / m
  Tmat = U_bar + (1 + 1/m) * B
  se = sqrt(diag(Tmat))
  
  df = (m - 1) * (1 + diag(U_bar) / ((1 + 1/m) * diag(B)))^2
  
  ci1 = Q_bar[1] + c(-1, 1) * qt(0.975, df[1]) * se[1]
  ci2 = Q_bar[2] + c(-1, 1) * qt(0.975, df[2]) * se[2]
  
  data.frame(
    est_CR = Q_bar[1], est_SR = Q_bar[2],
    bias_CR = Q_bar[1] - true1,
    bias_SR = Q_bar[2] - true2,
    cover_CR = as.numeric(ci1[1] <= true1 && true1 <= ci1[2]),
    cover_SR = as.numeric(ci2[1] <= true2 && true2 <= ci2[2]),
    pci_lb_CR = ci1[1], pci_ub_CR = ci1[2],
    pci_lb_SR = ci2[1], pci_ub_SR = ci2[2]
  )
}


#################################
# scenarios: full simulation ####
#################################

scenarios = function(condition_grid, minCR, maxCR, 
                     minSR, maxSR, lower, upper, iter = 100, 
                     true1 = 3, true2 = 3) {
  
  results_all = lapply(1:nrow(condition_grid), function(i) {
    cond = condition_grid[i, ]
    sim(
      seed = cond$seed,
      target = cond$target,
      meanCR = cond$meanCR,
      meanSR = cond$meanSR,
      sdCR = cond$sdCR,
      sdSR = cond$sdSR,
      minCR = minCR,
      maxCR = maxCR,
      minSR = minSR,
      maxSR = maxSR,
      lower = lower,
      upper = upper,
      iter = iter,
      true1 = true1,
      true2 = true2,
      distribution = cond$distribution
    )
    
  })
  do.call(rbind, results_all)
  
}

##################################
# sim: data generation + MNAR ####
##################################

sim = function(seed = NULL,
               target, meanCR, meanSR, sdCR, sdSR,
               minCR, maxCR, minSR, maxSR, lower, upper,
               iter = iter, true1 = 3, true2 = 3,
               distribution = c("normal", "uniform", "tmvn")) {
  
  distribution = match.arg(distribution)
  if (!is.null(seed)) set.seed(seed)
  
  S = 50
  N = 100
  Mu = c(0, 0)
  Tau = c(1, 1)
  rho = 0.7
  Sigma = diag(Tau) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(Tau)
  RTher = rmvnorm(S, Mu, Sigma)
  
  b0 = rnorm(S, 20, 2)
  b1 = rnorm(S, 0.5, 0.2)
  b2 = rnorm(S, 1.5, 0.3)
  b3 = 3
  
  data = lapply(1:S, function(i) {
    minA = runif(1, 18, 65)
    maxA = runif(1, minA + 5, 80)
    Age = runif(N, minA, maxA)
    Sex = rbinom(N, 1, 0.45)
    Ther = rbinom(N, 1, 0.5)
    
    CR = rnorm(N, b0[i] + b1[i] * Age + b2[i] * Sex + (b3 + RTher[i, 1]) * Ther, 5)
    SR = rnorm(N, b0[i] + b1[i] * Age + b2[i] * Sex + (b3 + RTher[i, 2]) * Ther, 7)
    
    data.frame(Study = i, Age, Sex = factor(Sex), Ther = factor(Ther), CR, SR)
  })
  
  d = do.call(rbind, data)
  
  dat = lapply(1:S, function(s) {
    Sn = d[d$Study == s, ]
    fitsur = systemfit(list(CR = CR ~ Age + Sex + Ther, SR = SR ~ Age + Sex + Ther), 
                       "SUR", data = Sn)
    sum = summary(fitsur)
    data.frame(
      Study = s,
      N = N,
      EstCR = sum$coefficients[4, 1],
      EstSR = sum$coefficients[8, 1],
      SECR = sum$coefficients[4, 2],
      SESR = sum$coefficients[8, 2],
      Cor.ws = sum$residCor["CR", "SR"]
    )
  })
  
  dat = do.call(rbind, dat)
  
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
  
  dat[M_cr == 1, c("EstCR", "SECR")] = NA
  dat[M_sr == 1, c("EstSR", "SESR")] = NA
  dat$Cor.ws[is.na(dat$EstCR) | is.na(dat$EstSR)] = NA
  
  res = genimp.mnar(
    df = dat,
    distribution = distribution,
    iter = iter,
    minCR = minCR, maxCR = maxCR,
    minSR = minSR, maxSR = maxSR,
    meanCR = meanCR, meanSR = meanSR,
    sdCR = sdCR, sdSR = sdSR,
    imprho = 0.7,
    meantmv = c(meanCR, meanSR),
    sigmatmv = matrix(c(sdCR^2, 
                        0.7 * sdCR * sdSR, 
                        0.7 * sdCR * sdSR, 
                        sdSR^2), nrow = 2),
    lower = lower,
    upper = upper
  )
  
  res_sum = sum.meth(res,
                      true1  = true1,
                      true2  = true2)   
  
  res_sum$replicate    = seed
  res_sum$distribution = distribution
  res_sum$meanCR       = meanCR
  res_sum$meanSR       = meanSR
  res_sum$target       = target
  
  return(res_sum)
}

######################
# conditions grid ####
######################

data_grid = expand.grid(
  seed = 1:5, # Irene: to be extended to 5000 or 10^4 these are the MC reps
  distribution = c("uniform", "normal", "tmvn"),
  target = c(0.10, 0.20, 0.30),
  stringsAsFactors = FALSE
)

uniform_tmvn = subset(data_grid, distribution != "normal")
uniform_tmvn$meanCR = ifelse(uniform_tmvn$distribution == "uniform", 0, 1)
# i set meanCR for multivariate normal to 1 to reflect that 'smaller' studies
# have been deleted

uniform_tmvn$sdCR = ifelse(uniform_tmvn$distribution == "uniform", 0, 6)
uniform_tmvn$meanSR = ifelse(uniform_tmvn$distribution == "uniform", 0, 5)
# i set meanCR for multivariate normal to 5 to reflect that 'bigger' studies
# have been deleted

uniform_tmvn$sdSR = ifelse(uniform_tmvn$distribution == "uniform", 0, 6)

normeans = c(0, 3, -3)
normal_grid = do.call(rbind, lapply(normeans, function(m) {
  g = subset(data_grid, distribution == "normal")
  g$meanCR = m
  g$meanSR = m
  g$sdCR = 6
  g$sdSR = 6
  g
}))

grid = rbind(uniform_tmvn, normal_grid)
head(grid)

##########################
# Run full simulation ####
##########################

library(parallel)

minCR = -100; maxCR = 100
minSR = -100; maxSR = 100
lower = c(-Inf, -Inf)
upper = c( Inf,  Inf)
iter  = 100
true1 = 3; true2 = 3


cl = makeCluster(detectCores() - 1)

clusterEvalQ(cl, {
  library(mixmeta); library(mvtnorm); library(systemfit); library(tmvtnorm)
})

clusterExport(
  cl,
  varlist = ls(),         
  envir   = .GlobalEnv
)

results_all_list <- parLapply(cl, 1:nrow(grid), function(i) {
  cond <- grid[i, ]
  sim(
    seed = cond$seed,
    target = cond$target,
    meanCR = cond$meanCR,
    meanSR = cond$meanSR,
    sdCR = cond$sdCR,
    sdSR = cond$sdSR,
    minCR = minCR,
    maxCR = maxCR,
    minSR = minSR,
    maxSR = maxSR,
    lower = lower,
    upper = upper,
    iter = iter,
    true1 = true1,
    true2 = true2,
    distribution = cond$distribution
  )
})

stopCluster(cl)

results_all = do.call(rbind, results_all_list)

results_all

################################
# Plot & summary of results ####
################################

summary_plot = results_all %>%
  group_by(distribution, meanCR, meanSR, target) %>%
  summarise(
    mean_estCR = mean(est_CR),
    lowerCR = mean(pci_lb_CR),
    upperCR = mean(pci_ub_CR),
    mean_estSR= mean(est_SR),
    lowerSR = mean(pci_lb_SR),
    upperSR = mean(pci_ub_SR),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(distribution, "(", meanCR, ",", meanSR, ")"))

ggplot(summary_plot, aes(x = mean_estCR, y = label)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lowerCR, xmax = upperCR), height = 0.25) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "#A4031F") +
  facet_wrap(~target) +
  labs(title = "Mean CI for CR across simulation replicates",
       x = "Estimated CR", y = "Imputation Method") +
  theme_minimal()

ggplot(summary_plot, aes(x = mean_estSR, y = label)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lowerSR, xmax = upperSR), height = 0.25) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "#A4031F") +
  facet_wrap(~target) +
  labs(title = "Mean CI for CR across simulation replicates",
       x = "Estimated CR", y = "Imputation Method") +
  theme_minimal()

results_summary = results_all %>%
  group_by(distribution, meanCR, meanSR, target) %>%
  summarise(
    mean_bias_CR = mean(bias_CR),
    mean_bias_SR = mean(bias_SR),
    coverage_CR = mean(cover_CR),
    coverage_SR = mean(cover_SR),
    .groups = "drop"
  ) 

results_summary = results_summary %>%
  mutate(fill_grp = case_when(
    distribution == "uniform" ~ "uniform",
    distribution == "tmvn"    ~ "mean = 1",      # for CR plot
    TRUE                      ~ paste0("mean = ", meanCR)
  ))


my_cols = c("mean = -3" = "#A2D729",
             "mean = 0"  = "#33A1FD",
             "mean = 3"  = "#840032",
             "uniform"   = "#F79824",
             "mean = 1"  = "#172A3A")   

ggplot(results_summary,
       aes(x = distribution,
           y = mean_bias_CR,
           fill = fill_grp)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(target ~ .) +
  scale_fill_manual(values = my_cols, name = NULL) +
  labs(title = "Mean bias of CR by distribution, mean and missingness percentage",
       x = NULL, y = "Mean bias (CR)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.spacing.y = unit(1, "lines"))


results_summary <- results_summary %>%
  mutate(fill_grp = case_when(
    distribution == "uniform" ~ "uniform",
    distribution == "tmvn"    ~ "mean = 5",      # for SR plot
    TRUE                      ~ paste0("mean = ", meanSR)
  ))

my_cols["mean = 5"] <- "#5C5198"

ggplot(results_summary,
       aes(x = distribution,
           y = mean_bias_SR,
           fill = fill_grp)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(target ~ .) +
  scale_fill_manual(values = my_cols, name = NULL) +
  labs(title = "Mean bias of SR by distribution, mean and missingness percentage",
       x = NULL, y = "Mean bias (SR)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.spacing.y = unit(1, "lines"))


# Plot coverage
ggplot(results_summary, aes(x= distribution, y = coverage_CR, fill = fill_grp)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = my_cols, name = NULL) +
  labs(title = "Coverage for CR", x = "Distribution", y = "Coverage") +
  facet_grid(~ target) +
  theme_minimal()

ggplot(results_summary, aes(x= distribution, y = coverage_SR, fill = fill_grp)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = my_cols, name = NULL) +
  labs(title = "Coverage for SR", x = "Distribution", y = "Coverage") +
  facet_grid(~ target) +
  theme_minimal()

# Irene: this can be refined, hic sunt leones
################################################################################
# Plot all panels
library(patchwork)

library(tidyr)

results_long = results_all[,c("replicate", "distribution", "target", "meanCR", "est_CR", "est_SR")]

avg_curve <- results_long %>%
  group_by(distribution, meanCR, target) %>%
  summarise(mean_est_CR = mean(est_CR), .groups = "drop")


plot_dat <- results_long %>% 
  mutate(fill_grp = ifelse(distribution == "uniform",
                           "uniform",
                           paste0("mean = ", meanCR)),
         target   = factor(target))     # make target discrete

ggplot(plot_dat,
       aes(x = target,
           y = est_CR,
           fill = fill_grp)) +
  geom_violin(trim = FALSE,
              scale = "width",
              colour = "grey30",
              alpha  = .5) +
  stat_summary(fun = mean, geom = "point",
               colour = "black", size = 2,
               position = position_dodge(.9)) +
  facet_wrap(~ distribution, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = my_cols,
                    name   = "mean") +
  labs(title = "Estimated CR",
       x = "Missingness percentage",
       y = "Estimated CR") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

ggplot(plot_dat,
       aes(x = target,
           y = est_SR,
           fill = fill_grp)) +
  geom_violin(trim = FALSE,
              scale = "width",
              colour = "grey30",
              alpha  = .5) +
  stat_summary(fun = mean, geom = "point",
               colour = "black", size = 2,
               position = position_dodge(.9)) +
  facet_wrap(~ distribution, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = my_cols,
                    name   = "mean") +
  labs(title = "Estimated SR",
       x = "Missingness percentage",
       y = "Estimated SR") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))


# Monte Carlo evaluation summary

results_mc_summary <- results_all %>%
  group_by(distribution, meanCR, target) %>%
  summarise(
    mean_bias_CR = mean(bias_CR),
    se_bias_CR   = sd(bias_CR) / sqrt(n()),
    
    mean_bias_SR = mean(bias_SR),
    se_bias_SR   = sd(bias_SR) / sqrt(n()),
    
    coverage_CR = mean(cover_CR),
    se_coverage_CR = sqrt(coverage_CR * (1 - coverage_CR) / n()),
    
    coverage_SR = mean(cover_SR),
    se_coverage_SR = sqrt(coverage_SR * (1 - coverage_SR) / n()),
    
    .groups = "drop"
  ) %>%
  mutate(fill_grp = ifelse(distribution == "uniform",
                           "uniform",
                           paste0("mean = ", meanCR)))
results_mc_summary
