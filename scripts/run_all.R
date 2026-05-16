# run_all.R

source("scripts/Figure1_histology_qPCR.R")
source("scripts/Figure2_microbiome.R")
source("scripts/Figure3_flow_cytometry.R")
source("scripts/Figure4_EAE.R")

sink("sessionInfo.txt")
sessionInfo()
sink()
