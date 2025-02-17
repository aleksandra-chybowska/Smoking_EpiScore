library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")
# datadir <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/data/"
# localdir <- "<local_data_dir>/methylation/GS_20k/Chromosomes/" # p17
# results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_cigs_day_w1_w3_w4/"
# phenotype <- "<cluster_home_dir>/Smoking/BayesRR/data/cigs_day17597_ASSIGN_W1_W3_W4.csv"
# residualise <- TRUE
# raw_data <- FALSE
# gs_wave <- "W1_W3_W4_cigs_day"
# pheno <- "cigs_day"
datadir <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/data/"
# results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_w1_w3_w4_complete_genetic"
results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final"
# phenotype <- "<cluster_home_dir>/Smoking/BayesRR/data/pack_years_17833_complete_genetic.csv"
phenotype <- "<cluster_home_dir>/Smoking/BayesRR/data/pack_years_17865_complete.csv"
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
savehistory(¨<cluster_home_dir>/Code/general/projects/smoking/summary_BayesR_pheno.R¨)
savehistory("<cluster_home_dir>/Code/general/projects/smoking/summary_BayesR_pheno.R")
