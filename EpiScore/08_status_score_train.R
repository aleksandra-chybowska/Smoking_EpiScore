library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")

datadir <- "<cluster_home_dir>/Smoking/Status_Score/data/"
results <- "<cluster_home_dir>/Smoking/Status_Score/results/unpaired/"
phenotype <- "<cluster_home_dir>/Smoking/BayesRR/data/pack_years_17865_complete.csv"
ever_smoke <- "<filespace_marioni_group_dir>/GS/GS_dataset/updated_smoking_jan_2019/ever_smoke.csv"
residualise <- TRUE
raw_data <- TRUE
gs_wave <- "W1_W3_W4"
pheno <- "ever_smoke"

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

## Import data
##########################################################################
if (raw_data == TRUE) {
  meth <- readRDS("<filespace_marioni_group_dir>/GS/GS_methylation/GS20k/mvals.rds") # p17
  target <- readRDS("<filespace_marioni_group_dir>/GS/GS_methylation/GS20k/GS20k_Targets.rds")

  meth <- subset(meth, colnames(meth) %in% target$Sample_Sentrix_ID)
  gc() # 831733  18413

  probes <- read.table("<filespace_marioni_group_dir>/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)$V1
  meth <- subset(meth, rownames(meth) %in% probes) # 752722  18869
  gc()

  # Pre-filtering by weight
  weights <- read.csv("<cluster_home_dir>/Smoking/Elnet_EpiScore/data/weights/Joehannes_pack_years_sup_tbl3.csv")
  print(dim(weights))
  weights <- weights[c("Name", "P.value")] # 18760     2 # Used to say "effect" here
  #weights <- weights[which(weights$P.value < 1e-4),] # 10414     2
  meth <- subset(meth, rownames(meth) %in% weights$Name) #  8982 18869
  gc()

  saveRDS(meth, "<cluster_home_dir>/Smoking/Status_Score/data/GS20k_mvals_filtered_by_weight.RDS", compress=F)

} else {
  target <- readRDS("<filespace_marioni_group_dir>/GS/GS_methylation/GS20k/GS20k_Targets.rds")
  meth <- readRDS("<cluster_home_dir>/Smoking/Status_Score/data/GS20k_mvals_filtered_by_weight.RDS")
}

ever_smoke <- read.csv(ever_smoke) # 24078
pheno <- merge(target, ever_smoke, by.x="Sample_Name", by.y="Sample_Name")
pheno = subset(pheno, ever_smoke != 5)

pdf("<cluster_home_dir>/Smoking/Status_Score/results/first/ever_smoke.pdf")
hist(pheno$ever_smoke)
dev.off()

pheno$status = ifelse(pheno$ever_smoke == 1, 1, ifelse(pheno$ever_smoke == 4, 3, 2))
rownames(pheno) = pheno$Sample_Sentrix_ID
write.csv(pheno, "<cluster_home_dir>/Smoking/Status_Score/data/ever_smoke_17730.csv")
pheno <- pheno[which(pheno$Sample_Sentrix_ID %in% colnames(meth)), ] # 17730

sample_size <- nrow(pheno)
print(paste0("Training in ", gs_wave, ": ", sample_size))

# Filter meth
meth <- meth[, pheno$Sample_Sentrix_ID] # 8982 17758
identical(pheno$Sample_Sentrix_ID, colnames(meth))
gc()

table(pheno$Sample_Sentrix_ID==colnames(meth))
table(rownames(pheno)==colnames(meth))
table(is.na(meth))


if(residualise == T) {
  design.resid <- model.matrix(~as.factor(sex) + age + as.factor(Set) , data=pheno)
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

saveRDS(meth, paste0(results, "meth_before_glmnet.RDS"), compress=F)
saveRDS(pheno, paste0(results, "pheno_before_glmnet.RDS"), compress=F)

## Elnet
##########################################################################
meth = readRDS(paste0(results, "meth_before_glmnet.RDS"))
pheno = readRDS(paste0(results, "pheno_before_glmnet.RDS"))

seed <- 42
folds <- 10

# x <- meth
y <- pheno$status

cv <- cv.glmnet(meth, y, family = "multinomial", type.multinomial = "ungrouped", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(meth, y, family = "multinomial", type.multinomial = "ungrouped", alpha = 0.5, lambda = cv$lambda.min)

temp <- coef(fit, s = fit$lambda) 
beta <- Reduce(cbind, temp)
beta <- beta[apply(beta != 0, 1, any),]
colnames(beta) <- names(temp)
beta <- as.matrix(beta)
beta <- as.data.frame(beta)
no_samples <- nrow(beta)

ret <- paste0(results, "/weights_", no_samples, "_samples.csv")
print(ret)
write.csv(beta, ret)
saveRDS(coef, ret)
saveRDS(cv, paste0(results,"/cv.RDS"))
saveRDS(fit, paste0(results,"/model.RDS"))


############# Doodles
