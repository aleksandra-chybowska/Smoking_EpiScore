library("ggplot2")
library("stringr")
library("pROC")
library("data.table")

epismoker <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/sensitivity/EpiSmokEr.csv")
results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final/"
pheno <- read.csv(paste0(results, "lbc36_predictions_pack_years_new_target.tsv"), sep = '\t')
pheno <- merge(pheno, epismoker[c("SampleName", "PredictedSmokingStatus")], by.x="basename", by.y="SampleName")
pheno$status_pred <- ifelse(pheno$PredictedSmokingStatus == "Current Smoker", 2, ifelse(pheno$PredictedSmokingStatus == "Never Smoker", 0, 1))

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
auc_36 <- auc(roc36_current_never)


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
