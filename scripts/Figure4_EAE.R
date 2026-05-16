# Figure4_EAE.R

source("scripts/00_setup.R")

# Input: data_processed/figure4_eae_scores.csv
# Required columns: mouse_id, group, day, clinical_score
eae <- read_csv("data_processed/figure4_eae_scores.csv") %>%
  mutate(group = factor(group, levels = c("Vehicle", "Mucinase")), day = as.numeric(day), mouse_id = factor(mouse_id))

eae_aov <- aov(clinical_score ~ group * day + Error(mouse_id / day), data = eae)
capture.output(summary(eae_aov), file = "outputs/statistics/Figure4B_EAE_two_way_ANOVA.txt")

eae_summary <- eae %>%
  group_by(group, day) %>%
  summarise(mean_score = mean(clinical_score, na.rm = TRUE), sem = sd(clinical_score, na.rm = TRUE) / sqrt(n()), n = n(), .groups = "drop")

write_csv(eae_summary, "outputs/statistics/Figure4B_EAE_clinical_score_summary.csv")

p_eae <- ggplot(eae_summary, aes(x = day, y = mean_score, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_score - sem, ymax = mean_score + sem), width = 0.3) +
  labs(x = "Days after immunization", y = "Clinical score", color = "Treatment", title = "EAE clinical disease course")

ggsave("outputs/figures/Figure4B_EAE_clinical_scores.pdf", p_eae, width = 5, height = 4)

peak_scores <- eae %>%
  group_by(mouse_id, group) %>%
  summarise(peak_score = max(clinical_score, na.rm = TRUE), .groups = "drop")

peak_stats <- t_test(peak_scores, peak_score ~ group, var.equal = FALSE) %>% add_significance("p")
write_csv(peak_stats, "outputs/statistics/Figure4B_peak_EAE_score_Welch_ttest.csv")

p_peak <- ggplot(peak_scores, aes(x = group, y = peak_score, color = group)) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  labs(x = NULL, y = "Peak clinical score", title = "Peak EAE clinical score") +
  theme(legend.position = "none")

ggsave("outputs/figures/Figure4B_peak_EAE_score.pdf", p_peak, width = 4, height = 4)

# Optional disease incidence: Kaplan-Meier/log-rank
# Input: data_processed/figure4_eae_incidence.csv
# Required columns: mouse_id, group, time_to_onset_or_censor, event
if (file.exists("data_processed/figure4_eae_incidence.csv")) {
  incidence <- read_csv("data_processed/figure4_eae_incidence.csv") %>%
    mutate(group = factor(group, levels = c("Vehicle", "Mucinase")), event = as.numeric(event))

  surv_obj <- Surv(time = incidence$time_to_onset_or_censor, event = incidence$event)
  fit <- survfit(surv_obj ~ group, data = incidence)
  logrank <- survdiff(surv_obj ~ group, data = incidence)

  capture.output(logrank, file = "outputs/statistics/Figure4B_EAE_incidence_logrank.txt")

  p_incidence <- ggsurvplot(
    fit,
    data = incidence,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Days after immunization",
    ylab = "Disease-free probability",
    legend.title = "Treatment"
  )

  ggsave("outputs/figures/Figure4B_EAE_incidence_KaplanMeier.pdf", p_incidence$plot, width = 5, height = 4)
}
