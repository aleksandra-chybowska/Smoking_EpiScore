#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena and Ola

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

datadir <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/data/"
localdir <- "<local_data_dir>/methylation/GS_20k/Chromosomes/" # p17
results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/joehannes/"

## Import data to test (W4)
##########################################################################

# Alcohol data
pack_years <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/pack_years_17865_complete.csv")
rownames(pack_years) <- pack_years$Sample_Sentrix_ID
pack_years <- pack_years[pack_years$Set == "wave4", ] # 8538

# Wave 4 methylation to test on
meth <- readRDS("<local_data_dir>/methylation/GS_20k/mvals.rds")
target <- readRDS("<local_data_dir>/methylation/GS_20k/GS20k_Targets_18869.rds")
meth <- meth[,target$X]
meth <- as.data.frame(meth)
gc()

# Filter meth
probes <- read.table("<filespace_marioni_group_dir>/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)$V1
anno <- readRDS("<filespace_marioni_group_dir>/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

meth <- meth[probes, pack_years$Sample_Sentrix_ID]
meth <- meth[rownames(meth) %in% rownames(common_anno),] # 386399 CpGs left
meth <- t(meth) #8538 386399
gc()

# Beta vals
meth <- m2beta(meth)


gc()


## Import predictors
##########################################################################

coefs_j <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/weights/Joehannes_categorical_sup_tbl2.csv")
coefs_j <- coefs_j[c("Probe.ID", "Effect")]

coefs_j_pack_years <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/weights/Joehannes_pack_years_sup_tbl3.csv")
coefs_j_pack_years <- coefs_j_pack_years[c("Name", "Effect")]

# cannot test here daniels probes as they were 
# coefs_t <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/weights/meanModelSM.csv")

rownames(coefs_j) <- coefs_j$Probe.ID #18760
rownames(coefs_j_pack_years) <- coefs_j_pack_years$Name #18760

#rownames(coefs_daniel) <- coefs_daniel$CpG
## Test in w4
##########################################################################

# Match CpGs
coefs_j <- coefs_j[which(rownames(coefs_j) %in% colnames(meth)),] # 16201/18760
coefs_j_pack_years <- coefs_j_pack_years[which(rownames(coefs_j_pack_years) %in% colnames(meth)),] # 16201/18760

meth_w4 <- meth[,rownames(coefs_j)]
meth_all_w4 <- meth[,rownames(coefs_j_pack_years)]

# Predictions
pred_w4 <- meth_w4 %*% coefs_j$Effect
pred_w4 <- as.data.frame(pred_w4)
names(pred_w4) <- c("py_pred")

pred_all_w4 <- meth_all_w4 %*% coefs_j_pack_years$Effect
pred_all_w4 <- as.data.frame(pred_all_w4)
names(pred_all_w4) <- c("py_pred")

# Export
write.table(data.frame(basename = rownames(pred_w4), pred_w4), paste0(results, "episcore_categorical.tsv"), sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_all_w4), pred_all_w4), paste0(results, "episcore_pack_years.tsv"), sep = "\t", row.names = F, quote = F)

# pred_w4  <- read.delim(paste0(results, "episcore_categorical.tsv"))
# pred_all_w4 <- read.delim(paste0(results, "episcore_pack_years.tsv"))

# Correlation?
pack_years <- pack_years[pred_w4$basename,]
pack_years$py_pred<- pred_w4$ac_pred
pack_years$py_pred_all <- pred_all_w4$ac_pred

r_w4 <- cor(pack_years$pack_years_clean, pack_years$py_pred, use="pairwise.complete.obs") # 0.619384 # trejo 0.700693700007292 #0.4945889
r_w4_all <- cor(pack_years$pack_years_clean, pack_years$py_pred_all, use="pairwise.complete.obs") # 0.619384 # trejo 0.700693700007292 #0.5019002
print(paste0("Correlation all waves: ", r_36))

# Incremental DNAm R2
null_36 <- summary(lm(pack_years_clean ~ age + sex, data=pheno))$r.squared # 0.03340952
full_36 <- summary(lm(pack_years_clean ~ age + sex + py_pred, data=pheno))$r.squared # 0.3960525
all_waves <- round(100*(full_36 - null_36), 3) # 36.264


r_w4 <- cor(alcohol$alcohol_units_log, alcohol$ac_pred, use="pairwise.complete.obs") # 0.394327
r_w4_all <- cor(alcohol$alcohol_units_log, alcohol$ac_pred_all, use="pairwise.complete.obs") # 0.4184333
r_w4_dan <- cor(alcohol$alcohol_units_log, alcohol$ac_pred_dan, use="pairwise.complete.obs") # 0.3346762

# Incremental R^2?
null_w4 <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol))$r.squared
full_w4 <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred, data=alcohol))$r.squared
full_w4_all <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_all, data=alcohol))$r.squared
full_w4_dan <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_dan, data=alcohol))$r.squared

round(100*(full_w4 - null_w4), 3) # 15.156
round(100*(full_w4_all - null_w4), 3) # 16.213
round(100*(full_w4_dan - null_w4), 3) # 10.685

# Incremental R^2? (Non-log)
null_w4 <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol))$r.squared
full_w4 <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred, data=alcohol))$r.squared
full_w4_all <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_all, data=alcohol))$r.squared
full_w4_dan <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_dan, data=alcohol))$r.squared

round(100*(full_w4 - null_w4), 3) # 15.671
round(100*(full_w4_all - null_w4), 3) # 16.17
round(100*(full_w4_dan - null_w4), 3) # 12.176


## Performance across categories (non-drinkers/moderate/heavy drinkers) AUC in LBC
###################################################################################

lbc_target <- read.table("<filespace_marioni_group_dir>/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T, row.names = 1)

# Binarize alcohol consumption, 0 = no drinkers/light drinkers, 1 = heavy drinkers
lbc_target$alcohol_cat_bi <- ifelse(lbc_target$alcohol_cat == "moderate-heavy_drinker", 1, 0)

# Import predictions
lbc_predictions_36 <- read.table("<filespace_marioni_group_dir>/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments.tsv", header = T)
lbc_predictions_21 <- read.table("<filespace_marioni_group_dir>/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments.tsv", header = T)

# Separate true values into 1936 and 1921
lbc_target_36 <- lbc_target[lbc_target$cohort == "LBC36",] # 2797   24
lbc_target_21 <- lbc_target[lbc_target$cohort == "LBC21",] # 692  24

# Match rownames
lbc_target_36 <- lbc_target_36[lbc_predictions_36$basename,]
lbc_target_21 <- lbc_target_21[lbc_predictions_21$basename,]

# ROC curves
roc_36 <- roc(response = lbc_target_36$alcohol_cat_bi, predictor = lbc_predictions_36$ac_pred)
roc_21 <- roc(response = lbc_target_21$alcohol_cat_bi, predictor = lbc_predictions_21$ac_pred)

auc_36 <- auc(roc_36) # 0.5942
auc_21 <- auc(roc_21) # 0.6531

## CHAT WITH DANIEL TO SEE HOW HE DID THIS!

## PLOTS