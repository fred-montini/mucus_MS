# 00_setup.R
# Setup script for mucus_MS analyses

required_packages <- c(
  "tidyverse",
  "ggplot2",
  "ggpubr",
  "vegan",
  "DESeq2",
  "survival",
  "survminer",
  "emmeans",
  "rstatix"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("DESeq2")
}

invisible(lapply(setdiff(required_packages, "DESeq2"), install_if_missing))

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(DESeq2)
library(survival)
library(survminer)
library(emmeans)
library(rstatix)

theme_set(theme_classic(base_size = 12))

dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/statistics", showWarnings = FALSE)
