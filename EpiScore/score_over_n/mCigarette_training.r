library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")

datadir <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/data/"
results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/score_over_n/"
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

  meth <- meth[,target$X]
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
  meth = readRDS('<cluster_home_dir>/Smoking/Elnet_EpiScore/data/20k_mvals_filtered_by_weight_1e-4tmp_full_epic_new_target.RDS')
}

pack_years <- read.csv(phenotype)
rownames(pack_years) <- pack_years$Sample_Sentrix_ID
pack_years <- pack_years[which(pack_years$Sample_Sentrix_ID %in% colnames(meth)), ] # 17758
# pack_years <- subset(pack_years, Set==gs_wave)
sample_size <- nrow(pack_years)

if (pheno == "cigs_day") {
  pack_years["py_clean_resid"] <- resid(lm(pack_years$cigs_day_clean ~ pack_years$age + factor(pack_years$sex), na.action = na.exclude)) 
} else {
  pack_years["py_clean_resid"] <- resid(lm(pack_years$pack_years_clean ~ pack_years$age + factor(pack_years$sex), na.action = na.exclude)) 
}
pack_years$py_clean_resid_scaled = scale(pack_years$py_clean_resid)

print(paste0("Training in ", gs_wave, ": ", sample_size))

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


iterations <- 1:17 * 1000
iterations <- c(iterations, nrow(meth))
sets_path <- paste0(datadir, "/training")
dir.create(sets_path)
big_meth = meth
big_pack_years = pack_years

saveRDS(big_meth, paste0(sets_path, "Big_meth.R"), compress = F)
saveRDS(big_pack_years, paste0(sets_path, "Big_pack_years.R"), compress = F)

for (it in iterations) {
  rows <- base::sample(rownames(big_meth), it, replace=F)
  meth = big_meth[rownames(big_meth) %in% rows,]
  pack_years = big_pack_years[rownames(big_pack_years) %in% rows,]
  rownames <- rownames(meth)
  
  identical(pack_years$Sample_Sentrix_ID, rownames(meth))  
  table(pack_years$Sample_Sentrix_ID==rownames(meth))
  table(rownames(pack_years)==rownames(meth))
  table(is.na(meth))

  #scale only at the end
  meth <- apply(meth, 2, scale) # 17758  9039
  rownames(meth) <- rownames
  gc()


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

  names(coefs)[1] <- "Coefficient" # Tidy naming  
  coefs$CpG_Site <- rownames(coefs) # Create cpg column
  coefs <- coefs[c(2,1)] # order, 1166  2 
  no_samples <- nrow(coefs)

  ret <- paste0(results, "weights_", it, "k.csv")
  print(ret)
  write.csv(coefs, ret, row.names = F)

}
