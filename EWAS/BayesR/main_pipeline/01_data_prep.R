#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena, modified by Ola

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")

############################################################################
# Script to prep files for BayesR+ (methylation)
############################################################################

datadir <- "<cluster_home_dir>/Smoking/BayesRR/data/"
localdir <- "<local_data_dir>/methylation/GS_20k/Chromosomes/" # p17
methdir <- "methylation_complete"
filename <- "pack_years_17865_complete.csv"
type <- "complete"

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

## (0) BMI file prep - sensitivity analysis, optional
check_sensitivity = FALSE

if (check_sensitivity == TRUE) {
  ds <- read.csv(paste0(datadir, "pack_years_bmi_17833.csv"))
  ds["bmi_resid"] <- resid(lm(ds$bmi ~ ds$age + factor(ds$sex), na.action = na.exclude)) 
  rownames(ds) = pack_years$Sample_Sentrix_ID

  write.table(x = t(as.matrix(as.numeric(scale(ds$bmi_resid)))),
  file = paste0(datadir, "bmi_17833_residualized.csvphen"), quote = F, sep = ",", row.names = F, col.names = F)

  saveRDS(ds, file = paste0(datadir, "ds_17833_residualized.rds"), compress=F)
}

## (1) Pheno file prep
############################################################################

# Import prepped pack_years data
pack_years <- read.csv(paste0(datadir, filename))
sample_size <- nrow(pack_years)
pack_years["pack_years_log_resid"] <- resid(lm(pack_years$pack_years_clean ~ 
  pack_years$age + factor(pack_years$sex), na.action = na.exclude)) 

rownames(pack_years) = pack_years$Sample_Sentrix_ID
# Export residualized
write.table(x = t(as.matrix(as.numeric(scale(pack_years$pack_years_log_resid)))),
file = paste0(datadir, "pack_years_", type, "_", sample_size, "_residualized.csvphen"),
              quote = F, sep = ",", row.names = F, col.names = F)

saveRDS(pack_years, file = paste0(datadir, "pack_years_", type, "_",
        sample_size, "_residualized.rds"), compress=F)

print(paste("Pheno file ready:", 
        paste0(datadir, "pack_years_", type, "_", sample_size, "_residualized.csvphen")))
# this is where I finished. I need to make Sextrix_IDs rownames and re-save.

## (2) Methylation files prep
############################################################################

# Set up stuff for paralellizing
cores <- detectCores() # 128
cl <- makeCluster(3, outfile=paste0(datadir, "parallel_test.txt")) # Let's start with 6 since they use up a lot of RAM
registerDoParallel(cl)

# Iterate per chromosome
probes <- read.table("<filespace_marioni_group_dir>/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)
foreach(i=1:22, .packages = "data.table") %dopar% { 
  print(paste0("Working on chromosome ",i))

  # Import
  meth <- readRDS(paste0(localdir, "GS20k_chr", i, "_mvals.rds"))
  
  # Subset
  meth <- meth[,which(colnames(meth) %in% rownames(pack_years))] # Subset to those with phenotype in question 
  meth <- meth[which(rownames(meth) %in% probes$V1),] # Subset to probes passing QC 
  
  # Match order of IDs in phenotype and methylation file 
  meth <- meth[,match(rownames(pack_years), colnames(meth))]
  
  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth, 1, meanimpute)
  meth <- t(meth)
  
  table(pack_years$Sample_Sentrix_ID==colnames(meth))
  table(rownames(pack_years)==colnames(meth))

  # Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~as.factor(sex) + age + as.factor(Batch), data=pack_years)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)),]
  rm(fit.resid)
  gc()
  
  # Write out CpGs 
  cpgs <- as.data.frame(row.names(meth))
  names(cpgs)[1] <- "CpG"
  fwrite(cpgs, paste0(datadir, methdir, "/GS20k_chr", i, "_cpgs.txt"),row.names=F)

  # Save out residualised file without scaling
  fwrite(meth, paste0(datadir, methdir, "/GS20k_chr", i, "_resid_mvals_notscaled.txt"),row.names=F)
  
  # Scale 
  meth <- t(apply(meth,1,scale))

  # Save out residualised file scaled
  fwrite(meth, paste0(datadir, methdir, "/GS20k_chr", i, "_resid_mvals.txt"),row.names=F)
  
  # Remove methylation object and clean up environment 
  rm(meth)
  gc()
} 

# End parallel
stopCluster(cl)

# Fuse methylation files - scaled
files <- list.files(path = paste0(datadir, methdir, "/"), pattern = "_resid_mvals.txt") # Get files
files <- files[order(files)]
data <- rbindlist(lapply(paste0(datadir, methdir, "/", files),fread))
gc()

# Export fused methylation file - scaled
fwrite(x = as.matrix(data), paste0(datadir, methdir, "/GS20k_allchrom.csv"), sep = ",", row.names = F, col.names = F, quote = F) 

# Fuse methylation files - not scaled
files <- list.files(path = paste0(datadir, methdir, "/"), pattern = "_resid_mvals_notscaled.txt") # Get files
files <- files[order(files)]
data <- rbindlist(lapply(paste0(datadir, methdir, "/", files),fread))
gc()

# Export fused methylation file - not scaled
fwrite(x = as.matrix(data), paste0(datadir, methdir, "/GS20k_allchrom_notscaled.csv"), sep = ",", row.names = F, col.names = F, quote = F) 

# Extract CpGs - will need this for processing final results files as the row.names get lost in methylation file in BayesR
cpgs <- list.files(path = paste0(datadir, methdir, "/"), pattern = "_cpgs.txt")
cpgs <- cpgs[order(cpgs)]

# Export CpGs
cg <- rbindlist(lapply(paste0(datadir, methdir, "/", cpgs),fread))
names(cg)[1] <- "Marker"
fwrite(x =as.matrix(cg), paste0(datadir, methdir, "/GS20k_allchrom_cpg_list.txt"), sep = ",", row.names = F, col.names = T, quote = F) 

## (4) Covariates file prep
############################################################################

# Daniel doesn't use any so maybe can get away without after adjusting for age/sex/batch/smoking episcore?
# "If you know about something, it's better to fit it as a covariate"

# Read in WBCs
wbc <- read.table(paste0(datadir, "wbc_quant.qcov"), header=T)
# Match order of IDs with other files 
wbc=wbc[match(pack_years$Sample_Sentrix_ID,wbc$FID),]
## Remove Sample ID information 
wbc$FID=NULL
wbc$IID=NULL 
wbc[1:ncol(wbc)]=apply(wbc[,1:ncol(wbc)],2,scale)

# Write out covariate file
write.table(x = as.matrix(wbc),
file = paste0(datadir, "basic_pack_years_cov_", type, ".csv"), quote = F, sep = ",", row.names = F, col.names = F)
