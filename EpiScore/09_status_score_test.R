library("ggplot2")
library("stringr")
library("pROC")
library("data.table")
setwd("<cluster_home_dir>/Smoking/Status_Score/")

results <- "results/first/"
run <- "1300_samples"
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

## Predictions in LBC
####################################################################

beta <- read.csv(paste0(results, "weights_", run, ".csv"))
rownames(beta) <- beta$X
beta = beta[-1]
beta <- beta[2:nrow(beta),] #1299    4

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
lbc_target_36 <- subset(lbc_target, cohort == "LBC36" & WAVE == 1) # 895  18
lbc_mvals_36 <- lbc_mvals[rownames(lbc_target_36),] # 895 1272

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

pred = NA
for (i in 1:ncol(beta)) {
# LBC36
	pred_36 <- lbc_mvals_36 %*% beta[,i]
	pred_36 <- as.data.frame(pred_36)
	names(pred_36) <- i

	pred = cbind(pred, pred_36)
}
pred[1] = rownames(pred)
colnames(pred)[1] = "Basename"

lbc_target_36 = merge(lbc_target_36, pred, by="Basename")

lbc_target_36$ID = as.character(lbc_target_36$ID)
pheno <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
pheno <- merge(lbc_target_36,
			pheno[c("lbc36no", "ever_smoke", "pack_years_clean")],
			by.x = "ID", by.y = "lbc36no")

pheno$status <- ifelse(pheno$ever_smoke == 2, 1, ifelse(pheno$ever_smoke == 0, 3, 2))
epismoker <- read.csv("../BayesR_EpiScore/data/sensitivity/EpiSmokEr.csv")
pheno <- merge(pheno, epismoker[c("SampleName", "PredictedSmokingStatus")], by.x="Basename", by.y="SampleName")
pheno$status_pred <- ifelse(pheno$PredictedSmokingStatus == "Current Smoker", 1, ifelse(pheno$PredictedSmokingStatus == "Never Smoker", 3, 2))

rownames(pheno) <- pheno$Basename
write.table(data.frame(basename = rownames(pheno), pheno), 
"predicted_status.tsv", sep = "\t", row.names = F, quote = F)


## Test performance
####################################################################
fit = readRDS("model.RDS")

pheno %>%
  mutate(dupa = predict(fit, lbc_mvals_36, type = "response")) %>%
  group_by(status) %>%
  summarize(class_rate = mean(pred == status))

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


epismoker <- read.csv("../BayesR_EpiScore/data/sensitivity/EpiSmokEr.csv")
results <- "results/first/"
pheno <- read.csv("data/lbc36_predictions_pack_years_new_target.tsv", sep="\t")
pheno <- merge(pheno, epismoker[c("SampleName", "PredictedSmokingStatus")], by.x="basename", by.y="SampleName")
pheno$status_pred <- ifelse(pheno$PredictedSmokingStatus == "Current Smoker", 1, ifelse(pheno$PredictedSmokingStatus == "Never Smoker", 3, 2))
pheno$status <- ifelse(pheno$ever_smoke == 2, 1, ifelse(pheno$ever_smoke == 0, 3, 2))

roc_to_df <- function(mod, roc) {
  ret = data.frame(
    "Model" = mod, 
    "Sensitivity" = roc$sensitivities,
    "Specificity" = roc$specificities
    )
  return(ret)
}

## Test performance
####################################################################

pheno_w1 <- pheno[pheno$WAVE == 1,] # 882
status_info <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
pheno_w1 <- merge(pheno_w1, status_info[c("lbc36no", "smokcat_w1")], by.x = "ID", by.y = "lbc36no")

current_never <- subset(pheno_w1, smokcat_w1 %in% c(0, 2))
former_never <- subset(pheno_w1, smokcat_w1 %in% c(0, 1))
current_former <- subset(pheno_w1, smokcat_w1 %in% c(1, 2))
# current - 2
# former - 1 
# never - 0

current_never$current <- ifelse(current_never$smokcat_w1 == 2, 1, 0)
former_never$former <- former_never$smokcat_w1
current_former$current <- ifelse(current_former$smokcat_w1 == 2, 1, 0)

roc36_current_never <- roc(response = current_never$current, predictor = current_never$status_pred)
auc_36 <- auc(roc36_current_never) # 0.9831


roc36_former_never <- roc(response = former_never$former, predictor = former_never$status_pred)
auc_36_former_never <- auc(roc36_former_never) # Area under the curve: 0.8499 (if tested with new target 0.8495)

roc36_current_former <- roc(response = current_former$current, predictor = current_former$status_pred)
auc_36_current_former <- auc(roc36_current_former) # Area under the curve: 0.8947 (new target: 0.8946)

### compare it with models containing age and sex 
current_never_model = glm(current ~ age + factor(sex) + status_pred, family="binomial", data=current_never)
roc_current_never_model <- roc(response = current_never$current, predictor = predict(current_never_model))

former_never_model = glm(former ~ age + factor(sex) + status_pred, family="binomial", data=former_never)
current_former_model = glm(current ~ age + factor(sex) + status_pred, family="binomial", data=current_former)

#roc <- roc(response = former_never$former, predictor = predict(full))

metrics_current_never <- roc_to_df("Current vs Never", roc36_current_never)
metrics_current_former <- roc_to_df("Current vs Former", roc36_current_former)
metrics_former_never <- roc_to_df("Former vs Never", roc36_former_never)
metrics <- rbind(metrics_current_never, metrics_current_former, metrics_former_never)

cbPalette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#F0E442" )

pdf(file=paste0(results, "AUC_Trejo.pdf"), height = 4, width = 6)

metrics %>%
  ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
  geom_line() +
  theme_light() +
  scale_colour_manual(values=cbPalette) +
  theme(legend.position = "right") +
  xlab("False Positive Rate (1 - Specificity)") +
  ylab("True Postive Rate (Sensitivity)")

dev.off()
