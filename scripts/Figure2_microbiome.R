# Figure2_microbiome.R

source("scripts/00_setup.R")

# Input: data_processed/figure2_species_relative_abundance.csv
# Required columns: sample_id, group, species, relative_abundance
species_abund <- read_csv("data_processed/figure2_species_relative_abundance.csv") %>%
  mutate(group = factor(group, levels = c("HC", "MS")), species = factor(species))

# Panel A: Shannon alpha diversity
species_wide <- species_abund %>%
  select(sample_id, species, relative_abundance) %>%
  pivot_wider(names_from = species, values_from = relative_abundance, values_fill = 0)

species_matrix <- species_wide %>% column_to_rownames("sample_id") %>% as.matrix()
shannon <- vegan::diversity(species_matrix, index = "shannon")

alpha_df <- tibble(sample_id = names(shannon), shannon_index = as.numeric(shannon)) %>%
  left_join(species_abund %>% distinct(sample_id, group), by = "sample_id") %>%
  mutate(group = factor(group, levels = c("HC", "MS")))

alpha_stats <- wilcox_test(alpha_df, shannon_index ~ group) %>% add_significance("p")
write_csv(alpha_stats, "outputs/statistics/Figure2A_Shannon_MannWhitney.csv")

p_alpha <- ggplot(alpha_df, aes(x = group, y = shannon_index, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  labs(x = NULL, y = "Shannon diversity index", title = "Fecal microbiota alpha diversity") +
  theme(legend.position = "none")

ggsave("outputs/figures/Figure2A_Shannon_alpha_diversity.pdf", p_alpha, width = 4, height = 4)

# Panel B: phylum-level relative abundance
# Input: data_processed/figure2_phylum_relative_abundance.csv
# Required columns: sample_id, group, phylum, relative_abundance
phylum_abund <- read_csv("data_processed/figure2_phylum_relative_abundance.csv") %>%
  mutate(group = factor(group, levels = c("HC", "MS")))

phylum_summary <- phylum_abund %>%
  group_by(group, phylum) %>%
  summarise(mean_relative_abundance = mean(relative_abundance, na.rm = TRUE), .groups = "drop")

write_csv(phylum_summary, "outputs/statistics/Figure2B_phylum_mean_relative_abundance.csv")

p_phylum <- ggplot(phylum_summary, aes(x = group, y = mean_relative_abundance, fill = phylum)) +
  geom_col(position = "fill") +
  labs(x = NULL, y = "Mean relative abundance", fill = "Phylum", title = "Average phylum-level composition")

ggsave("outputs/figures/Figure2B_phylum_relative_abundance.pdf", p_phylum, width = 5, height = 4)

# Panel C: targeted species-level comparisons by Mann-Whitney U test
target_species <- c(
  "Akkermansia muciniphila",
  "Akkermansia glycaniphila",
  "Bacteroides cellulosilyticus",
  "Phoenicola coprophilus",
  "Bifidobacterium breve",
  "Clostridium butyricum"
)

target_df <- species_abund %>% filter(species %in% target_species)

target_stats <- target_df %>%
  group_by(species) %>%
  wilcox_test(relative_abundance ~ group) %>%
  add_significance("p")

write_csv(target_stats, "outputs/statistics/Figure2C_target_species_MannWhitney.csv")

p_targets <- ggplot(target_df, aes(x = group, y = relative_abundance, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.8) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = NULL, y = "Relative abundance", title = "Targeted mucus-associated bacterial species") +
  theme(legend.position = "none", strip.text = element_text(face = "italic"))

ggsave("outputs/figures/Figure2C_target_species_relative_abundance.pdf", p_targets, width = 10, height = 6)

# Optional broad DESeq2 differential abundance analysis if integer count data are available
# Inputs: data_processed/microbiome_counts_species.csv and data_processed/microbiome_metadata.csv
if (file.exists("data_processed/microbiome_counts_species.csv") && file.exists("data_processed/microbiome_metadata.csv")) {
  counts_raw <- read_csv("data_processed/microbiome_counts_species.csv")
  metadata <- read_csv("data_processed/microbiome_metadata.csv")

  count_matrix <- counts_raw %>% column_to_rownames("taxon") %>% as.matrix()
  metadata <- metadata %>%
    filter(sample_id %in% colnames(count_matrix)) %>%
    arrange(match(sample_id, colnames(count_matrix))) %>%
    mutate(group = factor(group, levels = c("HC", "MS"))) %>%
    column_to_rownames("sample_id")

  count_matrix <- count_matrix[, rownames(metadata)]

  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = metadata, design = ~ group)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("group", "MS", "HC"))
  res_df <- as.data.frame(res) %>% rownames_to_column("taxon") %>% arrange(padj)
  write_csv(res_df, "outputs/statistics/Figure2_DESeq2_broad_taxonomic_differential_abundance.csv")
}
