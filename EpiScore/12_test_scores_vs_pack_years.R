library(dplyr)
# Check if missing and install if so 
list.of.packages <- c("data.table", "stringr","pROC", "ggplot2", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

# Load packages
library("ggplot2")
library("stringr")
library("pROC")
library("data.table")
# Set working directory
setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Replication")

#####################################################
######## 1. Create functions for tidying data #######
#####################################################
lbc_meth <- "/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/"
smoking_score <- read.csv(paste0(lbc_meth, "LBC_clock_output_3489.csv"))
smoking_min <- smoking_score[c("Basename", "Smoking")]

pheno <- read.csv("smoking_scores.csv")
ds <- merge(pheno, smoking_min, by = "Basename")
ds <- rename(ds, Grimage = Smoking)
colnames(ds)
write.csv(ds, "smoking_scores_grimage.csv", row.names = F)

####################################################################################
###### SCRIPT 2 : CALCULATE SCORE #5 - PCA+ELNET SCORE based on Wielscher CpGs  ####
####################################################################################

# Set working directory
#setwd("filepath/")

###############################################
####### 1. Load requisite libraries ###########
###############################################

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

# Read in previous scores
pheno = read.csv("smoking_scores_grimage.csv")
scores = c("mCigarette", "McCartney.et.al.",
           "Trejo.et.al.", "EpiSmokEr", "Grimage")

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

cat(report, file="continuous_scores_grimage.txt")

## ROC - only for binary outcomes
####################################################################

# pheno = na.omit(pheno)

# current_never <- subset(pheno, smokcat_w1 %in% c(0, 2))
# former_never <- subset(pheno, smokcat_w1 %in% c(0, 1))
# current_former <- subset(pheno, smokcat_w1 %in% c(1, 2))

# current_never$current <- ifelse(current_never$smokcat_w1 == 2, 1, 0)
# former_never$former <- former_never$smokcat_w1
# current_former$current <- ifelse(current_former$smokcat_w1 == 2, 1, 0)

# out = data.frame("score" = scores,  
#                  "current_never" = -1,
#                  "current_former" = -1,
#                  "former_never" = -1)

# for (score in scores) {
  
#   roc36_current_never <- roc(response = current_never$current, predictor = current_never[,score])
#   auc_36 <- auc(roc36_current_never)
#   out[which(out$score == score), "current_never"] = auc_36
  
#   # current_never_model = glm(current ~ mCigarette, family="binomial", data=current_never)
#   # roc_current_never_model <- roc(response = current_never$current, predictor = predict(current_never_model))
#   # auc <- auc(roc_current_never_model)
#   # auc
  
#   roc36_former_never <- roc(response = former_never$former, predictor = former_never[,score])
#   auc_36_former_never <- auc(roc36_former_never) # Area under the curve: 0.8499 (if tested with new target 0.8495)
#   out[which(out$score == score), "former_never"] = auc_36_former_never
  
#   roc36_current_former <- roc(response = current_former$current, predictor = current_former[,score])
#   auc_36_current_former <- auc(roc36_current_former) # Area under the curve: 0.8947 (new target: 0.8946)
#   out[which(out$score == score), "current_former"] = auc_36_current_former
# }

# write.table(out, "AUC_table_grimage.csv", sep=',')

