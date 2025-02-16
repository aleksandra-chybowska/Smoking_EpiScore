#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("caret")

############################################################################
# Script to prep files for BayesR+ (genetic)
############################################################################

datadir <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/"
methdir <- "methylation_W3_W1_no_age_sex"
gendir <- "genetics_W1_W3"
filename <- "pack_years_9316_wave1_wave3_genetic_age_sex.csv"
type <- "wave1_wave3_genetic_age_sex"

meanimpute <- function(x) ifelse(is.na(x), mean(x,na.rm=T), x)


## (1) Phenotype file prep
############################################################################

# No need to re-prep, still have from methylation analysis.
# Daniel also residualizes with genetics PCs, but will skip this time around.


## (2) Genetics files prep
############################################################################

# Import alcohol file to make sure order is the same as pheno file
pack_years <- read.csv(paste0(datadir, filename))

gen <- as.data.frame(fread("/Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS.fam"))
pack_years <- pack_years[pack_years$Sample_Name %in% gen$V2,] # 16689

# Filter fam to just samples with alcohol data + methylation data
fam <- read.table("/GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015.fam")
fam <- fam[which(fam$V2 %in% pack_years$Sample_Name),]

# Export samples to filter genotype files with PLINK
fam_file = paste0(datadir, "pack_years_", type, "_fam.txt")
write.table(fam[,1:2], file = fam_file, sep = "\t", quote = F, row.names = F, col.names = F)
print("Prepped for PLINK!")
print(fam_file)

# In terminal 
# cd /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/genetics_W1_W3/chr
# for chr in {1..22}
# do
# plink19 --bfile /GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015 --chr $chr --make-bed --out pack_years_W1_W3_chr${chr}
# plink19 --recodeA --bfile pack_years_W1_W3_chr${chr} --keep ../../pack_yearswave1_wave3_genetic_age_sex_fam.txt --maf 0.01 --out pack_years_W1_W3_chr${chr}
# done

# Set up stuff for paralellizing
cores <- detectCores() # 128
cl <- makeCluster(3, outfile=paste0(datadir,"parallel_test.txt")) # Let's start with 6 since they use up a lot of RAM
registerDoParallel(cl)

# Prep genotype data, per chromosome and in parallel
foreach(i=1:22, .packages = "data.table") %dopar% { 
    print(paste0("Working on chromosome ",i))

    # Import genotype data and stratify ahead of imputation/scaling
    snps <- fread(paste0(datadir, gendir, "/chr/pack_years_W1_W3_chr", i,".raw"), header=T, stringsAsFactors=F)
    snps <- as.data.frame(snps)
    rownames(snps) <- as.numeric(snps$IID)
    snps <- snps[,7:ncol(snps)]
	ids  <- match(pack_years$Sample_Name, rownames(snps))
    snps <- snps[ids,] # Match order to pheno
	identical(as.integer(rownames(snps)), pack_years$Sample_Name)


    # Mean impute
    snps <- apply(snps, 2, meanimpute)
    gc()

    # Scale
    snps <- apply(snps, 2, scale)
    gc()

    # Remove SNPs with NAs (0 variance SNPs)
    snps <- snps[,colSums(is.na(snps)) < nrow(snps)]
    gc()

    # Save SNP names
    snp_names <- as.data.frame(colnames(snps))
    names(snp_names)[1] <- "Marker"
    fwrite(x =as.matrix(snp_names), paste0(datadir, gendir, "/pack_years_snp_list_chr", i, ".txt"), 
	sep = ",", row.names = F, col.names = T, quote = F) 

    # Export
    snps <- t(snps)
    fwrite(x = as.matrix(snps), paste0(datadir, gendir, "/chr/pack_years_processed_chr", i, ".csv"), 
	sep = ",", row.names = F, col.names = F, quote = F) 

    rm(snps)
    gc()
}

# End parallel
stopCluster(cl)

# Fuse chromosomes
files <- list.files(path = paste0(datadir, gendir, "/chr/"), pattern = "pack_years_processed_chr") # Get files
files <- files[order(files)]
data <- rbindlist(lapply(paste0(datadir, gendir, "/chr/", files), fread))
fwrite(x = as.matrix(data), paste0(datadir, gendir, "pack_years_W3_W1_allchrom_processed.csv"), sep = ",", row.names = F, col.names = F, quote = F) 
rm(data)
gc()


# Fuse SNPs
snps <- list.files(path = paste0(datadir, gendir), pattern = "pack_years_snp_list_chr")
snps <- snps[order(snps)]
snps_all <- rbindlist(lapply(paste0(datadir, gendir, "/", snps), fread))
names(snps_all)[1] <- "Marker"
fwrite(x =as.matrix(snps_all), paste0(datadir, gendir, "/pack_years_snp_list.txt"), sep = ",", row.names = F, col.names = T, quote = F) 


## (3) Prep genetics + methylation file
############################################################################

# Sex-agnostic
snps <- fread(paste0(datadir, gendir, "pack_years_W3_W1_allchrom_processed.csv"))
snps <- as.data.frame(snps)
cpgs <- fread(paste0(datadir, "/methylation_W3_W1_no_age_sex/GS20k_allchrom.csv"))
cpgs <- as.data.frame(cpgs)

# add Age and Sex from pheno
agesex <- data.frame(age = log(pack_years$age), sex = pack_years$sex)
rownames(agesex) <- agesex$Sample_Name
agesex$sex <- ifelse(agesex$sex == "F", 1, 0)
agesex <- scale(agesex)
names <- as.data.frame(colnames(agesex))
colnames(names) <- c("Marker")

agesex <- t(agesex)
data <- rbind(snps, cpgs)
data <- rbind(data, agesex)

fwrite(x = data, paste0(datadir, "multiomics_W1_W3/pack_years_W3_W4_allchrom_processed.csv"), 
	sep = ",", row.names = F, col.names = F, quote = F) 
rm(snps, cpgs, data)
gc()

# Fuse SNPs and CpGs
snps_list <- read.csv(paste0(datadir, gendir, "/pack_years_snp_list.txt"))
cpgs_list <- read.csv(paste0(datadir, "methylation_W3_W1_no_age_sex/GS20k_allchrom_cpg_list.txt"))
markers <- rbind(snps_list, cpgs_list)
markers <- rbind(markers, names)

fwrite(x = as.matrix(markers), paste0(datadir, "multiomics_W1_W3/pack_years_snp_cpg_agesex_list.txt"), 
	sep = ",", row.names = F, col.names = T, quote = F) 


## (4) Prep genetics + methylation groups file
############################################################################

markers["group"] <- c(rep(0, length(snps_list$Marker)), 
					  rep(1, length(cpgs_list$Marker)),
					  rep(1, length(names$Marker)))
fwrite(x = as.matrix(markers), paste0(datadir, "multiomics_W1_W3/pack_years_cpg_groups.txt"), 
					  sep=' ', row.names = F, col.names = F, quote = F)