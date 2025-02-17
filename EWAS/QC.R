# load libraries
library(tidyverse)
library(QCEWAS)
library(data.table)

setwd("<cluster_home_dir>/Smoking/EWAS")

ont = read.csv("<cluster_home_dir>/ONT/EWAS/5_40_46_full_EWAS_ONT.csv") # need to update with full pvals
path = '<cluster_home_dir>/TWIST/Manual_EWAS/'

twist = read.csv(paste0(path, "pass10_40_46_full_EWAS_TWIST.csv"))
epic = read.csv(paste0(path, 'W3_46_processed_new_annots.csv'))

format = read.table("sample1.txt.gz")
# 1    PROBEID              BETA               SE                P_VAL
# 2 cg21517506 -26.6469042651649 8.76670790587349  0.00236926733677786
# 3 cg12080369   25.452315181525 7.08179149299185 0.000325576880336428

formatted_epic = epic[c("cpg", "beta", "se", "pval")]
colnames(formatted_epic) = c("PROBEID", "BETA", "SE", "P_VAL")

formatted_twist = twist[c("cpg", "beta", "se", "p")]
colnames(formatted_twist) = c("PROBEID", "BETA", "SE", "P_VAL")

formatted_ont = ont[c("cpg", "beta", "se", "p")]
colnames(formatted_ont) = c("PROBEID", "BETA", "SE", "P_VAL")

QC_results <- EWAS_QC(data = formatted_epic,
outputname = "formatted_epic_46",
save_final_dataset = FALSE)

QC_results_twist <- EWAS_QC(data = formatted_twist,
outputname = "formatted_twist_46",
save_final_dataset = FALSE)

QC_results_ont <- EWAS_QC(data = formatted_ont,
outputname = "formatted_ont_46",
save_final_dataset = FALSE)
