# load the data
res <- readRDS("sim0907.rds")
library(tidyr)

combined_df <- res %>%
  mutate(sim_id = row_number()) %>%
  unnest(res, names_sep = "_") %>%
  group_by(distribution, target, meanCR)

split_dfs <- combined_df %>%
  group_by(distribution, target, meanCR) %>%
  group_split()

restot <- bind_rows(split_dfs)
restot$method <- paste(restot$distribution, restot$meanCR, restot$meanSR, sep = ", ")

sum_plot = restot %>%    
  group_by(distribution, meanCR, meanSR, target, method) %>%
  summarise(
    mean_estCR = mean(res_est_CR),
    lowerCR    = mean(res_pci_lb_CR),
    upperCR    = mean(res_pci_ub_CR),
    mean_estSR = mean(res_est_SR),
    lowerSR    = mean(res_pci_lb_SR),
    upperSR    = mean(res_pci_ub_SR),
    .groups    = "drop"
  ) %>%
  mutate(label = paste0(distribution, "(", meanCR, "," , meanSR, ")"))

library(ggplot2)

ggplot(sum_plot, aes(x = factor(meanCR), y = mean_estCR, color = distribution, group = distribution)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(aes(ymin = lowerCR, ymax = upperCR), 
                position = position_dodge(width = 0.4), width = 0.2) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray40") +
  facet_wrap(~ target) +
  labs(
    title = "Estimated CR with 95% CI by Imputation Distribution and Missingness Percentage",
    x = "Mean of Imputation Distribution (meanCR)",
    y = "Estimated CR",
    color = "Distribution"
  ) +
  theme_minimal()

ggplot(sum_plot, aes(x = factor(meanSR), y = mean_estSR, color = distribution, group = distribution)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(aes(ymin = lowerSR, ymax = upperSR), 
                position = position_dodge(width = 0.4), width = 0.2) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray40") +
  facet_wrap(~ target) +
  labs(
    title = "Estimated SR with 95% CI by Imputation Distribution and Missingness Percentage",
    x = "Mean of Imputation Distribution (meanSR)",
    y = "Estimated SR",
    color = "Distribution"
  ) +
  theme_minimal()

