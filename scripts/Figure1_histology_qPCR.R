# Figure1_histology.R

source("scripts/00_setup.R")

# Input: data_processed/figure1_histology.csv
# Required columns: sample_id, group, measure, value
histology <- read_csv("data_processed/figure1_histology.csv") %>%
  mutate(group = factor(group, levels = c("HC", "MS")), measure = factor(measure))

histology_stats <- histology %>%
  group_by(measure) %>%
  t_test(value ~ group, var.equal = TRUE) %>%
  add_significance("p")

write_csv(histology_stats, "outputs/statistics/Figure1_histology_ttests.csv")

p_histology <- ggplot(histology, aes(x = group, y = value, color = group)) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  facet_wrap(~ measure, scales = "free_y") +
  labs(x = NULL, y = "Quantification", title = "Histological quantification") +
  theme(legend.position = "none")

ggsave("outputs/figures/Figure1B_histology_quantification.pdf", p_histology, width = 8, height = 4)

