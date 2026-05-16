# Figure3_flow_cytometry.R

source("scripts/00_setup.R")

# Input: data_processed/figure3_flowcytometry_frequencies.csv
# Required columns: sample_id, group, panel, population, frequency
flow <- read_csv("data_processed/figure3_flowcytometry_frequencies.csv") %>%
  mutate(group = factor(group, levels = c("HC", "MS")), panel = factor(panel), population = factor(population))

flow_stats <- flow %>%
  group_by(panel, population) %>%
  wilcox_test(frequency ~ group) %>%
  add_significance("p")

write_csv(flow_stats, "outputs/statistics/Figure3_flowcytometry_MannWhitney.csv")

p_flow <- ggplot(flow, aes(x = group, y = frequency, color = group)) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  facet_wrap(panel ~ population, scales = "free_y") +
  labs(x = NULL, y = "Frequency (%)", title = "Intestinal immune cell subset frequencies") +
  theme(legend.position = "none", strip.text = element_text(size = 8))

ggsave("outputs/figures/Figure3_flowcytometry_frequencies.pdf", p_flow, width = 12, height = 8)
