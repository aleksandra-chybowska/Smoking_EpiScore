library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")
# datadir <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/data/"
# localdir <- "/Local_Data/methylation/GS_20k/Chromosomes/" # p17
# results <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/results/j_1e-4_cigs_day_w1_w3_w4/"
# phenotype <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/cigs_day17597_ASSIGN_W1_W3_W4.csv"
# residualise <- TRUE
# raw_data <- FALSE
# gs_wave <- "W1_W3_W4_cigs_day"
# pheno <- "cigs_day"
datadir <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/data/"
# results <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_w1_w3_w4_complete_genetic"
results <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final"
# phenotype <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_17833_complete_genetic.csv"
phenotype <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_17865_complete.csv"
residualise <- TRUE
raw_data <- TRUE
gs_wave <- "W1_W3_W4"
pheno <- "pack_years"
pack_years <- read.csv(phenotype)
summary(pack_years)
hist(pheno$age)
hist(pack_years$age)
sd(pack_years$age)
pack_years$sex = as.factor(pack_years$sex)
summary(pack_years)
prop.table(pack_years$sex)
prop.table(table(pack_years$sex))
10551/(10551+7314)
savehistory(¨/Cluster_Filespace/Marioni_Group/Ola/Code/general/projects/smoking/summary_BayesR_pheno.R¨)
savehistory("/Cluster_Filespace/Marioni_Group/Ola/Code/general/projects/smoking/summary_BayesR_pheno.R")
