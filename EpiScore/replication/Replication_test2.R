####################################################################################
###### SCRIPT 2 : CALCULATE SCORE #5 - PCA+ELNET SCORE based on Wielscher CpGs  ####
####################################################################################

# Set working directory
#setwd("filepath/")

###############################################
####### 1. Load requisite libraries ###########
###############################################

# Check if missing and install if so 
list.of.packages <- c("data.table", "stringr","pROC", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library("ggplot2")
library("stringr")
library("pROC")
library("data.table")

roc_to_df <- function(mod, roc) {
  ret = data.frame(
    "Model" = mod, 
    "Sensitivity" = roc$sensitivities,
    "Specificity" = roc$specificities
  )
  
  return(ret)
}

#################################################################
#### 1. Read in scores from Script 1 - to merge in at end  ######
#################################################################

setwd("<cluster_home_dir>/Smoking/Replication")
# Read in previous scores
pheno = read.csv("smoking_scores.csv")
scores = c("mCigarette", "McCartney.et.al.", "Trejo.et.al.", "EpiSmokEr")

## Test performance
####################################################################
report = ""

for (score in scores) {

  report = paste0(report, sprintf("Score: %s \n", score))
  
  r <- cor(pheno$pack_years_clean, pheno[score], use="pairwise.complete.obs")
  report = paste0(report, sprintf("Correlation with pack years: %.3f \n", r))
  
  # Incremental DNAm R2
  null <- summary(lm(pheno$pack_years_clean ~ pheno$age + pheno$sex))$r.squared
  full <- summary(lm(pheno$pack_years_clean ~ pheno$age + pheno$sex + pheno[,score]))$r.squared 
  w1 <- round(100*(full - null), 3) # 41.861
  
  report = paste0(report, sprintf("Null R2: %.3f \n", null))
  report = paste0(report, sprintf("Full R2: %.3f \n", full))
  report = paste0(report, sprintf("Incremental R2: %.3f \n", w1))
}

cat(report, file="continuous_scores.txt")

## ROC
####################################################################

pheno = na.omit(pheno)

current_never <- subset(pheno, smokcat_w1 %in% c(0, 2))
former_never <- subset(pheno, smokcat_w1 %in% c(0, 1))
current_former <- subset(pheno, smokcat_w1 %in% c(1, 2))

current_never$current <- ifelse(current_never$smokcat_w1 == 2, 1, 0)
former_never$former <- former_never$smokcat_w1
current_former$current <- ifelse(current_former$smokcat_w1 == 2, 1, 0)

out = data.frame("score" = scores,  
                 "current_never" = -1,
                 "current_former" = -1,
                 "former_never" = -1)

for (score in scores) {
  # score = "mCigarette"
  roc36_current_never <- roc(response = current_never$current, predictor = current_never[,score])
  auc_36 <- auc(roc36_current_never)
  out[which(out$score == score), "current_never"] = auc_36
  
  # current_never_model = glm(current ~ mCigarette, family="binomial", data=current_never)
  # roc_current_never_model <- roc(response = current_never$current, predictor = predict(current_never_model))
  # auc <- auc(roc_current_never_model)
  # auc
  
  roc36_former_never <- roc(response = former_never$former, predictor = former_never[,score])
  auc_36_former_never <- auc(roc36_former_never) # Area under the curve: 0.8499 (if tested with new target 0.8495)
  out[which(out$score == score), "former_never"] = auc_36_former_never
  
  roc36_current_former <- roc(response = current_former$current, predictor = current_former[,score])
  auc_36_current_former <- auc(roc36_current_former) # Area under the curve: 0.8947 (new target: 0.8946)
  out[which(out$score == score), "current_former"] = auc_36_current_former

  # metrics_current_never <- roc_to_df("Current vs Never", roc36_current_never)
  # metrics_current_former <- roc_to_df("Current vs Former", roc36_current_former)
  # metrics_former_never <- roc_to_df("Former vs Never", roc36_former_never)
  # metrics <- rbind(metrics_current_never, metrics_current_former, metrics_former_never)
  # write.table(metrics, "Figure_3_source_data_mCigarette_metrics.csv", sep=',')

}

write.table(out, "AUC_table.csv", sep=',')

# #Basic boxplot
# library(tidyverse)
# library(hrbrthemes)
# library(viridis)

# pdf("boxplots_for_Riccardo.pdf")
# pheno %>%
#   ggplot( aes(x=smokcat_w1, y=mCigarette, fill=as.factor(smokcat_w1))) +
#     geom_boxplot() +
#     ggtitle("Basic boxplot") +
#     xlab("")
# dev.off()


#roc <- roc(response = former_never$former, predictor = predict(full))

# metrics_current_never <- roc_to_df("Current vs Never", roc36_current_never)
# metrics_current_former <- roc_to_df("Current vs Former", roc36_current_former)
# metrics_former_never <- roc_to_df("Former vs Never", roc36_former_never)
# metrics <- rbind(metrics_current_never, metrics_current_former, metrics_former_never)
# 
# cbPalette <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#F0E442" )
# 
# pdf(file=paste0(results, "AUC_Trejo.pdf"), height = 4, width = 6)
# 
# metrics %>%
#   ggplot( aes(x=1-Specificity, y=Sensitivity, group=Model, color=Model)) +
#   geom_line() +
#   theme_light() +
#   scale_colour_manual(values=cbPalette) +
#   theme(legend.position = "right") +
#   xlab("False Positive Rate (1 - Specificity)") +
#   ylab("True Postive Rate (Sensitivity)")
# 
# dev.off()
