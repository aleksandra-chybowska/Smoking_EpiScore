#######################################
###### SCRIPT 2 : TEST SCORES 1-4  ####
#######################################

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
#### 1. Read in scores from Script 1  ######
#################################################################
# setwd("<path_to_replication_folder>")
# Read in previous scores
pheno <- read.csv("results/smoking_scores.csv")
scores <- c("mCigarette", "McCartney.et.al.", "Trejo.et.al.", "EpiSmokEr")

########################
## Test performance #### 
######################## 
report <- ""

for (score in scores) {

  report <- paste0(report, sprintf("Score: %s \n", score))
  
  r <- cor(pheno$pack_years, pheno[score], use="pairwise.complete.obs")
  report <- paste0(report, sprintf("Correlation with pack years: %.3f \n", r))
  
  # Incremental DNAm R2
  null <- summary(lm(pheno$pack_years ~ pheno$age + pheno$sex))$r.squared
  full <- summary(lm(pheno$pack_years ~ pheno$age + pheno$sex + pheno[,score]))$r.squared 
  w1 <- round(100*(full - null), 3) # 41.861
  
  report <- paste0(report, sprintf("Null R2: %.3f \n", null))
  report <- paste0(report, sprintf("Full R2: %.3f \n", full))
  report <- paste0(report, sprintf("Incremental R2: %.3f \n", w1))
}

cat(report, file="results/continuous_scores.txt")

###############
####  ROC  #### 
###############

pheno = na.omit(pheno)

# divide by smoking category: 2 - current smokers, 1 - former smokers, 0 - never smokers
current_never <- subset(pheno, smoking_cat %in% c(0, 2))
former_never <- subset(pheno, smoking_cat %in% c(0, 1))
current_former <- subset(pheno, smoking_cat %in% c(1, 2))

# divide categories into binary classes 
current_never$current <- ifelse(current_never$smoking_cat == 2, 1, 0)
former_never$former <- former_never$smoking_cat
current_former$current <- ifelse(current_former$smoking_cat == 2, 1, 0)

out <- data.frame("score" = scores,  
                 "current_never" = -1,
                 "current_former" = -1,
                 "former_never" = -1)

for (score in scores) {
  
  roc_current_never <- roc(response = current_never$current, predictor = current_never[,score])
  auc <- auc(roc_current_never)
  out[which(out$score == score), "current_never"] = auc
  
  roc_former_never <- roc(response = former_never$former, predictor = former_never[,score])
  auc_former_never <- auc(roc_former_never) # Area under the curve: 0.8499 (if tested with new target 0.8495)
  out[which(out$score == score), "former_never"] = auc_former_never
  
  roc_current_former <- roc(response = current_former$current, predictor = current_former[,score])
  auc_current_former <- auc(roc_current_former) # Area under the curve: 0.8947 (new target: 0.8946)
  out[which(out$score == score), "current_former"] = auc_current_former
}

write.csv(out, "results/AUC_data.csv", row.names = F)

