##################################################################################
### Script to generate DNAm predictor of Smoking in LBC1936 at wave 1 (age 70) ###
##################################################################################

### Look at phenotypes - choose outcome ###
library("foreign")
library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("caret")

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/data/")
run <- "wave3"
bayesR <- paste0("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/results/runs/", run, "/")
input <- paste0(bayesR, "summary/pack_years_4380_only_W3_meanbeta_pip.tsv")
output <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/runs/wave3/"

### Finished prepping pheno

### Generate methylation score ###

### Import methylation data ###

ph <- read.csv("LBC/pheno_min_1072.csv")
file <- input
beta <- read.table(file, row.names = 2, sep = "\t", header = TRUE)

# Import LBC data
lbc_target <- readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/targets_3489_bloodonly.rds")
rownames(lbc_target) <- lbc_target$Basename

lbc_mvals <- readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
lbc_mvals <- t(lbc_mvals)

# Match CpGs
beta <- beta[which(rownames(beta) %in% colnames(lbc_mvals)),]
lbc_mvals <- lbc_mvals[,rownames(beta)]

# Divide LBC data into LBC21 and LBC36
lbc_target_36 <- lbc_target[lbc_target["cohort"] == "LBC36",] # 2797, waves 1, 2, 3, and 4
lbc_mvals_36 <- lbc_mvals[rownames(lbc_target_36),]

# Scale M-vals
lbc_mvals_36 <- scale(lbc_mvals_36)

# Replace missing values with 0
lbc_mvals_36[is.na(lbc_mvals_36)] <- 0

## Predictions
#################################################################

# LBC36
pred_36 <- lbc_mvals_36 %*% beta$Mean_Beta
pred_36 <- as.data.frame(pred_36)

names(pred_36) <- c("py_pred")
lbc_target_36["py_pred"] <- pred_36["py_pred"]
pheno <- merge(lbc_target_36,
            ph[c("lbc36no", "ever_smoke", "pack_years_clean")],
            by.x = "ID", by.y = "lbc36no")
rownames(pheno) = pheno$Basename
write.table(data.frame(basename = rownames(pheno), pheno), paste0(output, 
"lbc36_predictions_pack_years_", run, ".tsv"), sep = "\t", row.names = F, quote = F)


## Assess performance of predictor 
#################################################################

# Correlation
r_36 <- cor(pheno$pack_years_clean, pheno$py_pred, use="pairwise.complete.obs") # 0.619384
print(paste0("Correlation all waves: ", r_36))
# Incremental DNAm R2
null_36 <- summary(lm(pack_years_clean ~ age + sex, data=pheno))$r.squared # 0.03340952
full_36 <- summary(lm(pack_years_clean ~ age + sex + py_pred, data=pheno))$r.squared # 0.3960525
all_waves <- round(100*(full_36 - null_36), 3) # 36.264

print(paste0("Null_36: ", null_36))
print(paste0("Full_36: ", full_36))
print(paste0("Incremental R2: ", all_waves))
# P-vals
summary(lm(pack_years_clean ~ age + sex, data = pheno))
summary(lm(pack_years_clean ~ age + sex + py_pred, data = pheno))

# Difference in alcohol consumption between top and bottom decile of prediction
pheno$py_pred_decile <- as.numeric(cut(pheno$py_pred, quantile(pheno$py_pred, prob = seq(0, 1, length = 11)), include.lowest=T ))
summary(lm(pack_years_clean ~ age + sex + factor(py_pred_decile), data = pheno))

# Fully adjusted incremental DNAm R2
# nullf <- summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung, data=w1_36))$r.squared
# fullf <- summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + g_pred, data=w1_36))$r.squared
# round(100*(fullf - nullf), 3)

# P-vals for full model
# summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + g_pred, data=w1_36))
# summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + X1.000000, data=w1_36))


## Assess performance of predictor - just wave 1
#################################################################

pheno_w1 <- pheno[pheno$WAVE == 1,] # 882

r_36_w1 <- cor(pheno_w1$pack_years_clean, pheno_w1$py_pred, use="pairwise.complete.obs") # 0.6560837
print(paste0("Correlation wave1: ", r_36_w1 ))
# Incremental DNAm R2
null_36_w1 <- summary(lm(pack_years_clean ~ age + sex, data=pheno_w1))$r.squared # 0.03018441
full_36_w1 <- summary(lm(pack_years_clean ~ age + sex + py_pred, data=pheno_w1))$r.squared # 0.4487904
w1 <- round(100*(full_36_w1 - null_36_w1), 3) # 41.861
print(paste0("Null_36_w1: ", null_36_w1))
print(paste0("Full_36_w1: ", full_36_w1))
print(paste0("Incremental R2 W1: ", w1))

# P-vals
summary(lm(pack_years_clean ~ age + sex, data = pheno_w1))
summary(lm(pack_years_clean ~ age + sex + py_pred, data = pheno_w1))

# Difference in alcohol consumption between top and bottom decile of prediction
pheno_w1$py_pred_decile <- as.numeric(cut(pheno_w1$py_pred, quantile(pheno_w1$py_pred, prob = seq(0, 1, length = 11)), include.lowest=T ))
summary(lm(pack_years_clean ~ age + sex + factor(py_pred_decile), data = pheno_w1))

