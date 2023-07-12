###############################################
#### BAYESR+ EWAS of Pack Years ###############
###############################################

## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/") 

## Load requisite libraries
library(data.table)
library(limma)

## Create function for mean imputation
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/"

#########################################################
######## STEP 1 - PREPARATION OF PHENOTYPE FILE ########
#########################################################

# Read in phenotype from OSCA (already prepared)
phen=read.table(paste0(path, "Phenotypes/pack_years.phen"), header=F)

# Read in with methylation basenames  
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

## Subset to those with complete genetic data 
gen=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS.fam"))
samps=merge(samps,gen,by.x="Sample_Name",by.y="V2")

# Merge phenotype data in with remaining basenames 
phenos=merge(phen,samps[,c("Sample_Name","Sample_Sentrix_ID")],by.x="V1",by.y="Sample_Name")

## Make additional IDs file where we have FID and IID 
fam=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")
fam1=merge(fam,phenos,by.x="volid", by.y="V1")
gwas.ids <- data.frame(FID =fam1$famid,
                       IID = fam1$volid)
write.table(gwas.ids,paste0(path, "BayesR/IDs/pack_years_ids.list"),row.names=F,quote=F,sep="\t")

# Write out phenotype 
write.table(x = t(as.matrix(as.numeric(scale(phenos[,3])))),
	file=paste0(path, "BayesR/Phenotypes/pack_years.csvphen"),quote = F, sep = ",", row.names = F, col.names = F)


##########################################################
######## STEP 2 - PREPARATION OF METHYLATION FILE ########
##########################################################

# Read in methylation file # by chromosome 
for(i in 1:22){ 
  meth=readRDS(paste0("GS20k_chr", i, "_mvals.rds"))
  
  # Subset to those with phenotype in question 
  meth=meth[,which(colnames(meth) %in% phenos$Sample_Sentrix_ID)]
  
  # Subset to probes passing QC 
  probes=read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)
  meth=meth[which(row.names(meth) %in% probes$V1),]
  
  # Match order of IDs in phenotype and methylation file 
  ids=phenos$Sample_Sentrix_ID
  meth=meth[,match(ids,colnames(meth))]
  
  # Check order of IDs match between phenotype and methylation files 
  table(colnames(meth)==phenos$Sample_Sentrix_ID)
  
  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth,1,meanimpute)
  
  # Transpose back original format - as apply changes the format
  meth=t(meth)
  
  # Prepare covariate matrix for regressions 
  
  # Match order of IDs with other files 
  samps=samps[match(ids,samps$Sample_Sentrix_ID),]
  
  # Check order of IDs match with other files 
  table(samps$Sample_Sentrix_ID==phenos$Sample_Sentrix_ID)
  table(samps$Sample_Sentrix_ID==colnames(meth))
  
  ## Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~sex + age + as.factor(Batch), data=samps)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)), ]
  rm(fit.resid)
  gc()
  
  ## Scale 
  meth=t(apply(meth,1,scale))
  
  ## Write out CpGs 
  cpgs=as.data.frame(row.names(meth))
  names(cpgs)[1]="CpG"
  fwrite(cpgs, paste0(path, "BayesR/CpGs/pack_years/GS20k_chr", i, "_cpgs.txt"),row.names=F)
  
  # Save out residualised file 
  fwrite(meth, paste0(path, "BayesR/Methylation/pack_years/GS20k_chr", i, "_resid_mvals.txt"),row.names=F)
  
  ## Remove methylation object and clean up environment 
  rm(meth)
  gc()
  print(i)
} 

## Write out ID order for later
saveRDS(ids, paste0(path, "BayesR/IDs/pack_years_meth_ids.rds"))
