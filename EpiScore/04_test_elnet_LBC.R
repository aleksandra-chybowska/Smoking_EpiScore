library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")

run <- "W1_W3_W4_1256_methylation_adjusted_set_correct_target"
results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final/"

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

## Predictions in LBC
####################################################################
sink(paste0(results, "predictor_", run, ".csv"))

beta <- read.csv(paste0(results, "weights_", run, ".csv"))
rownames(beta) <- beta$CpG_Site
beta <- beta[2:nrow(beta),] #1165    3

#intercept <- coefs[1,2]
#intercept_all <- coefs_all[1,2]

# Import LBC data
lbc_target <- readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/targets_3489_bloodonly.rds")
rownames(lbc_target) <- lbc_target$Basename
lbc_mvals <- readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
lbc_mvals <- t(lbc_mvals)
# lbc_mvals <- m2beta(lbc_mvals)

# Match CpGs
beta <- beta[which(rownames(beta) %in% colnames(lbc_mvals)),]
lbc_mvals <- lbc_mvals[,rownames(beta)]

# Divide LBC data into LBC21 and LBC36
lbc_target_36 <- lbc_target[lbc_target["cohort"] == "LBC36",] # 2797, waves 1, 2, 3, and 4
lbc_mvals_36 <- lbc_mvals[rownames(lbc_target_36),] # 2797 1145

# Scale M-vals
lbc_mvals_36 <- scale(lbc_mvals_36)

# Replace missing values with 0
lbc_mvals_36[is.na(lbc_mvals_36)] <- 0

# > identical(colnames(lbc_mvals_36), rownames(beta))
# [1] TRUE
# > identical(rownames(lbc_mvals_36), rownames(lbc_target_36))
# [1] TRUE

# Keep just wave 1
#lbc_target_36 <- lbc_target_36[lbc_target_36$WAVE == 1,] # 895


## Predictions
#################################################################
setwd(results)

# LBC36
pred_36 <- lbc_mvals_36 %*% beta$Coefficient
pred_36 <- as.data.frame(pred_36)

names(pred_36) <- c("py_pred")
lbc_target_36["py_pred"] <- pred_36["py_pred"]
lbc_target_36$ID = as.character(lbc_target_36$ID)
pheno <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
pheno <- merge(lbc_target_36,
            pheno[c("lbc36no", "ever_smoke", "pack_years_clean")],
            by.x = "ID", by.y = "lbc36no")

rownames(pheno) <- pheno$Basename
write.table(data.frame(basename = rownames(pheno), pheno), 
"lbc36_predictions_pack_years_new_target_elnet.tsv", sep = "\t", row.names = F, quote = F)

## Test performance
####################################################################

pheno_w1 <- pheno[pheno$WAVE == 1,] # 882

r_36_w1 <- cor(pheno_w1$pack_years_clean, pheno_w1$py_pred, use="pairwise.complete.obs") # 0.6560837
print(paste0("Correlation wave1: ", r_36_w1 ))
# Incremental DNAm R2
null_36_w1 <- summary(lm(pack_years_clean ~ age + sex, data=pheno_w1))$r.squared 
full_36_w1 <- summary(lm(pack_years_clean ~ age + sex + py_pred, data=pheno_w1))$r.squared 
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
sink()
