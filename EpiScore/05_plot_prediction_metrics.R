library("ggplot2")
library("stringr")
library("pROC")
library("data.table")
library("ggtext")

results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final/"
pheno <- read.csv(paste0(results, "lbc36_predictions_pack_years_new_target.tsv"), sep = '\t')
#pheno <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/runs/sensitivity/lbc36_predictions_pack_years_sensitivity_Daniel.tsv", sep="\t")
#pheno <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/runs/sensitivity/lbc36_predictions_pack_years_sensitivity_Trejo.tsv", sep="\t")
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

## ROC
####################################################################

status_info <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
pheno_w1 <- merge(pheno_w1, status_info[c("lbc36no", "smokcat_w1")], by.x = "ID", by.y = "lbc36no")

pdf(paste0(results, "cor_pack_years_py_episcore_Trejo.pdf"), 
width = 5, height = 4)
#plot(twist["chr5-373240-373241"], df[,"cg05575921"])
ggplot(pheno_w1, aes(x=pack_years_clean, y=py_pred)) +
  geom_point() +
  xlab("Pack years") +
  ylab("EpiScore") +
  geom_smooth(method='lm', colour = "red")
dev.off()
# measured pack years vs predicted pack years 

current_never <- subset(pheno_w1, smokcat_w1 %in% c(0, 2))
former_never <- subset(pheno_w1, smokcat_w1 %in% c(0, 1))
current_former <- subset(pheno_w1, smokcat_w1 %in% c(1, 2))

current_never$current <- ifelse(current_never$smokcat_w1 == 2, 1, 0)
former_never$former <- former_never$smokcat_w1
current_former$current <- ifelse(current_former$smokcat_w1 == 2, 1, 0)

roc36_current_never <- roc(response = current_never$current, predictor = current_never$py_pred)
auc_36 <- auc(roc36_current_never) # 0.9831
proc <- MLmetrics::PRAUC(y_pred = current_never$py_pred, y_true = current_never$current)

roc36_former_never <- roc(response = former_never$former, predictor = former_never$py_pred)
auc_36_former_never <- auc(roc36_former_never) # Area under the curve: 0.8499 (if tested with new target 0.8495)
proc_former_never <- MLmetrics::PRAUC(y_pred = former_never$py_pred, y_true = former_never$former)

roc36_current_former <- roc(response = current_former$current, predictor = current_former$py_pred)
auc_36_current_former <- auc(roc36_current_former) # Area under the curve: 0.8947 (new target: 0.8946)
proc_current_former <- MLmetrics::PRAUC(y_pred = current_former$py_pred, y_true = current_former$current)

### compare it with models containing age and sex 
current_never_model = glm(current ~ age + factor(sex) + py_pred, family="binomial", data=current_never)
roc_current_never_model <- roc(response = current_never$current, predictor = predict(current_never_model))

former_never_model = glm(former ~ age + factor(sex) + py_pred, family="binomial", data=former_never)
current_former_model = glm(current ~ age + factor(sex) + py_pred, family="binomial", data=current_former)

#roc <- roc(response = former_never$former, predictor = predict(full))

metrics_current_never <- roc_to_df("Current vs Never", roc36_current_never)
metrics_current_former <- roc_to_df("Current vs Former", roc36_current_former)
metrics_former_never <- roc_to_df("Former vs Never", roc36_former_never)
metrics <- rbind(metrics_current_never, metrics_current_former, metrics_former_never)

cbPalette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#F0E442" )

pdf(file=paste0(results, "AUC_Test.pdf"), height = 5, width = 5)

metrics %>%
  ggplot(aes(x = 1 - Specificity, y = Sensitivity, group = Model, color = Model)) +
  geom_line() +
  theme_light() +
  scale_color_manual(
      values = cbPalette,
      labels = c(paste("<span style='color:",  cbPalette[2],"'>Current vs Never: 0.98</span>"), 
      paste("<span style='color:",  cbPalette[1],"'>Current vs Former: 0.90</span>"), 
      paste("<span style='color:",  cbPalette[3],"'>Former vs Never: 0.85</span>")
    )
  ) +
  guides(color = guide_legend(
    override.aes = list(linetype = "blank", size = 0) # Removes legend lines
  )) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = c(0.75, 0.2), # Adjust position as needed
    legend.background = element_rect(fill = "white", color = "#afa393"),
    legend.text = element_markdown(color=cbPalette, size = 12),
    legend.key.width = unit(0, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.title = element_text(face = "bold") # Optional: bold title
  ) +
  labs(color = "AUC") + # Sets the legend title
  xlab("False Positive Rate (1 - Specificity)") +
  ylab("True Positive Rate (Sensitivity)")

dev.off()


# metrics %>%
#   ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
#   geom_line() +
#   theme_light() +
#   scale_colour_manual(values=cbPalette) +
#   theme(legend.position = "right") +
#   xlab("False Positive Rate (1 - Specificity)") +
#   ylab("True Postive Rate (Sensitivity)")