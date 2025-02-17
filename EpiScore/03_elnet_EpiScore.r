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
results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final_test"
# phenotype <- "<cluster_home_dir>/Smoking/BayesRR/data/pack_years_17833_complete_genetic.csv"
phenotype <- "<cluster_home_dir>/Smoking/BayesRR/data/pack_years_17865_complete.csv"
residualise <- TRUE
raw_data <- FALSE
gs_wave <- "W1_W3_W4"
pheno <- "pack_years"


m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

## Import data
##########################################################################
if (raw_data == TRUE) {
  meth <- readRDS("<local_data_dir>/methylation/GS_20k/mvals.rds") # p17
  #target <- readRDS("<local_data_dir>/methylation/GS_20k/GS20k_Targets_18869.rds")
  target <- readRDS("<filespace_marioni_group_dir>/GS/GS_methylation/GS20k/GS20k_Targets.rds")

  meth <- meth[,target$Sample_Sentrix_ID]
  meth <- as.data.frame(meth)
  gc() # 831733  18413

  #colnames - people 

  probes <- read.table("<filespace_marioni_group_dir>/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)$V1
  meth <- meth[probes, ] # 752722  18869

  # anno <- readRDS("<filespace_marioni_group_dir>/Daniel/EPIC_AnnotationObject_df.rds")
  # common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]
  # meth <- meth[rownames(meth) %in% rownames(common_anno), ] # 17758  9039 CpGs left
  gc()

  # Pre-filtering by weight
  weights <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/weights/Joehannes_pack_years_sup_tbl3.csv")
  print(dim(weights))
  weights <- weights[c("Name", "P.value")] # 18760     2 # Used to say "effect" here
  #weights <- weights[which(weights$P.value < 1e-4),] # 10414     2
  meth <- meth[rownames(meth) %in% weights$Name, ] #  8982 18869
  gc()

  saveRDS(meth, "<cluster_home_dir>/Smoking/Elnet_EpiScore/data/20k_mvals_filtered_by_weight_1e-4tmp_full_epic_new_target.RDS", compress=F)

} else {
  meth <- readRDS("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/20k_mvals_filtered_by_weight_1e-4tmp_full_epic_new_target.RDS")
  weights <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/weights/Joehannes_pack_years_sup_tbl3.csv")
  weights <- weights[c("Name", "Effect")] # 18760     2
  weights <- weights[which(weights$Effect < 1e-4),] # 10414     2
  meth <- meth[rownames(meth) %in% weights$Name, ] #  8982 18869
  meth <- na.omit(meth)
  gc()
}

pack_years <- read.csv(phenotype)
rownames(pack_years) <- pack_years$Sample_Sentrix_ID
pack_years <- pack_years[which(pack_years$Sample_Sentrix_ID %in% colnames(meth)), ] # 17758
# pack_years <- subset(pack_years, Set==gs_wave)
sample_size <- nrow(pack_years)
print(paste0("Training in ", gs_wave, ": ", sample_size))
if (pheno == "cigs_day") {
  pack_years["py_clean_resid"] <- resid(lm(pack_years$cigs_day_clean ~ pack_years$age + factor(pack_years$sex), na.action = na.exclude)) 
} else {
  pack_years["py_clean_resid"] <- resid(lm(pack_years$pack_years_clean ~ pack_years$age + factor(pack_years$sex), na.action = na.exclude)) 
}
pack_years$py_clean_resid_scaled = scale(pack_years$py_clean_resid)

# Filter meth
meth <- meth[, pack_years$Sample_Sentrix_ID] # 8982 17758
identical(pack_years$Sample_Sentrix_ID, colnames(meth))
gc()

table(pack_years$Sample_Sentrix_ID==colnames(meth))
table(rownames(pack_years)==colnames(meth))
table(is.na(meth))

# saveRDS(meth, paste0(datadir, "20k_before_res.RDS"), compress=F)
#meth <- readRDS("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/20k_before_res.RDSFALSE")

if(residualise == T) {
  design.resid <- model.matrix(~as.factor(sex) + age + as.factor(Set) , data=pack_years)
  fit.resid <- limma::lmFit(meth, design.resid) 
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)),]
  rm(fit.resid)
  gc() # 10350 17758
}

meth <- t(meth)
meth <- m2beta(meth)
table(is.na(meth))
rownames <- rownames(meth)
meth <- apply(meth, 2, scale) # 17758  9039
rownames(meth) <- rownames
gc()

saveRDS(meth, paste0(results, '/prep_meth.RDS'), compress = F)
saveRDS(pack_years, paste0(results, '/pheno.RDS'), compress = F)
## Elnet
##########################################################################
seed <- 42
folds <- 10

# x <- meth
y <- pack_years$py_clean_resid_scaled

cv <- cv.glmnet(meth, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(meth, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)

coefs <- coef(fit) # Extract coeficients  
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 

print(cv$lambda.min) #  Elnet: 0.012577, Lasso: 0.006288498, Ridge: 6.288498
print(dim(coefs)) # Elnet: 1256, Lasso, 1246, Ridge: 16202
names(coefs)[1] <- "Coefficient" # Tidy naming  
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order, 1166  2 
no_samples <- nrow(coefs)

ret <- paste0(results, "/weights_", gs_wave, "_", no_samples, "_methylation_adjusted_set_correct_target.csv")
print(ret)
write.csv(coefs, ret, row.names = F)

