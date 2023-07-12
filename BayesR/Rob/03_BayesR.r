
###############################################################
######## STEP 3 - PREPARE GENETIC FILES  ######################
###############################################################

####### Recode step - not in R ######

plink19 \
--bfile /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS \
--recode A \
--out /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/Genetics_Raw/recode \


### Back in R #####
library(data.table)

path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/"
phen=read.table(paste0(path, "Phenotypes/pack_years.phen"), header=F)

# Read in with methylation basenames  
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

## Subset to those with complete genetic data 
gen=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS.fam"))
samps=merge(samps,gen,by.x="Sample_Name",by.y="V2")

# Merge phenotype data in with remaining basenames 
phenos=merge(phen,samps[,c("Sample_Name","Sample_Sentrix_ID")],by.x="V1",by.y="Sample_Name")

## Read in resulting file in R for preparation step 
snp=as.data.frame(fread(paste0(path, "BayesR/Genetics_Raw/recode.raw")))

## Match order with other files
snp=snp[which(snp$IID %in% phenos$V2),]
ids1=phenos$V2
snp=snp[match(ids1,snp$IID),]

# Remove redundant columns
snp[,1:6]=NULL

# Extract SNP names 
snps=as.data.frame(colnames(snp))
names(snps)[1]="Marker"
fwrite(x =as.matrix(snps), paste0(path, "BayesR/SNPs/pack_years_snp_list.csv"), sep = ",", row.names = F, col.names = T, quote = F) 

# Mean impute any missing values 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
snp <- apply(snp,2,meanimpute)

# Transpose back original format - as apply changes the format
# snp=t(snp)

# Scale each value 
snp <- apply(snp,2,scale)

snp=t(snp)

# Save out file 
fwrite(x =as.matrix(snp), paste0(path, "BayesR/Inputs/Genetics/pack_years_17833.csv"), sep = ",", row.names = F, col.names = F, quote = F) 



###############################################################
######## STEP 4 - PREPARE FIXED EFFECT COVARIATE FILES ########
###############################################################

# Read in WBCs
wbc=read.table("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/wbc_quant.qcov", header=T)

# Combine with sample ID information - to help merge in with other covariates 
wbc1=merge(wbc,samps[,c("Sample_Sentrix_ID", "Sample_Name")],by.x="FID",by.y="Sample_Sentrix_ID")

# Combine with BMI 
bod=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/body.csv")
wbc1=merge(wbc1, bod[,c("id", "bmi")], by.x="Sample_Name", by.y="id")

# Mean impute BMI 
wbc1$bmi=meanimpute(wbc1$bmi)

# Log transform and scale BMI 
wbc1$bmi=log(wbc1$bmi)

# Match order of IDs with other files 
wbc1=wbc1[match(ids,wbc1$FID),]

saveRDS(wbc1, paste0(path, "BayesR/covars_imperfect.rds"), compress = F)
## Remove Sample ID information 
wbc1$FID=NULL
wbc1$IID=NULL 
wbc1$Sample_Name=NULL 
wbc1[1:ncol(wbc1)]=apply(wbc1[,1:ncol(wbc1)],2,scale)

# Write out covariate file
write.table(x = as.matrix(wbc1),file = paste0(path, "BayesR/Covariates/basic_pack_years.csv") ,quote = F, sep = ",", row.names = F, col.names = F)



##############################################################
######## STEP 5 - COMBINE INDIVIDUAL CHROMOSOME FILES ########
##############################################################

# Change working directory 
setwd(paste0(path, "BayesR/Methylation/pack_years/"))

# Extract files 
files=list.files(".", ".txt")
files=files[order(files)]

# Read in and rbind all methylation files 
data <- rbindlist(lapply(files,fread))

## Write out final file
fwrite(x =as.matrix(data), paste0(path, "BayesR/Inputs/Methylation/pack_years_17833.csv"), sep = ",", row.names = F, col.names = F, quote = F) 

#======= HERE ========

# Extract CpGs - will need this for processing final results files as the row.names get lost in methylation file in BayesR
cpgs=list.files("../../CpGs/pack_years/", ".txt")
cpgs=cpgs[order(cpgs)]

# Ensure that order is same as methylation file just created 
# methylation
ids.cpg=gsub("GS20k_chr", "", files)
ids.cpg=gsub("_.*", "", ids.cpg)
# cpgs
ids1.cpg=gsub("GS20k_chr", "", cpgs)
ids1.cpg=gsub("_.*", "", ids1.cpg)
# is order same? 
table(ids.cpg==ids1.cpg)

# Read in cpg lists 
# Change working directory 
setwd(paste0(path, "BayesR/CpGs/pack_years/"))
cg <- rbindlist(lapply(cpgs,fread))
names(cg)[1]="Marker"
# Save out file
fwrite(x =as.matrix(cg), paste0(path, "BayesR/CpGs/pack_years_cpg_list.csv"), sep = ",", row.names = F, col.names = T, quote = F) 



##############################################################
######## STEP 6 - COMBINE GENETICS + METHYLATION FILE ########
##############################################################

## Combine with genetics file to make multi-omics file 
new.data=rbind(data,snp)

# Remove other files to clear up environment
rm(data)
rm(snp)
gc()

## Write out final file
fwrite(x =as.matrix(new.data), paste0(path, "BayesR/Inputs/Combined/pack_years_17833.csv"), sep = ",", row.names = F, col.names = F, quote = F) 

# Remove object and clear up space
rm(new.data)
gc()

# Combine markers to make groups file 
cg$Group=0
snps$Group=1
new.markers=rbind(cg,snps)
# Save out file 
fwrite(x =as.matrix(new.markers), paste0(path, "BayesR/Groups/pack_years.txt"),  sep=' ', row.names = F, col.names = T, quote = F) 


#####################################
#### STEP 7 - BAYESR  #####
#####################################

## Epigenetics Alone 
cd /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/

../../../BayesRRcmd/src/brr --data-file Inputs/Methylation/pack_years_17833.csv --pheno Phenotypes/pack_years.csvphen --analysis-type preprocess --fixed_effects Covariates/basic_pack_years.csv --fixedEffectNumber 6 --thread 24 --thread-spawned 24 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Inputs/Methylation/pack_years_17833.csv --pheno Phenotypes/pack_years.csvphen --fixed_effects Covariates/basic_pack_years.csv --fixedEffectNumber 6 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.01,0.1,1.0" --mcmc-samples Outputs/Methylation/pack_years.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1

## Genetics Alone 

cd /Cluster_Filespace/Marioni_Group/Rob/Smoking/BayesR/
  
../../../BayesRRcmd/src/brr --data-file Inputs/Genetics/pack_years_17833.csv --pheno Phenotypes/pack_years.csvphen --analysis-type preprocess --fixed_effects Covariates/basic_pack_years.csv --fixedEffectNumber 6 --thread 24 --thread-spawned 24 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Inputs/Genetics/pack_years_17833.csv --pheno Phenotypes/pack_years.csvphen --fixed_effects Covariates/basic_pack_years.csv --fixedEffectNumber 6 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.0001,0.001,0.01" --mcmc-samples Outputs/Genetics/pack_years.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1


## Combined 

cd /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/
  
../../../BayesRRcmd/src/brr --data-file Inputs/Combined/pack_years_17833.csv \
--pheno Phenotypes/pack_years.csvphen \
--analysis-type preprocess \
--fixed_effects Covariates/basic_pack_years.csv \
--fixedEffectNumber 6 \
--group Groups/pack_years.txt \
--thread 24 \
--thread-spawned 24 --marker-cache --seed 1 > Outputs/Combined/combined.out


# ../../../BayesRRcmd/src/brr --data-file Inputs/Combined/pack_years_17833.csv --pheno Phenotypes/pack_years.csvphen --analysis-type preprocess --fixed_effects Covariates/basic_pack_years.csv --fixedEffectNumber 6 --group Groups/pack_years.txt --marker-cache --seed 1 > Outputs/Combined/combined.out

../../../BayesRRcmd/src/brr --data-file Inputs/Combined/pack_years_17833.csv --pheno Phenotypes/pack_years.csvphen --fixed_effects Covariates/basic_pack_years.csv --fixedEffectNumber 6 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --group Groups/pack_years.txt --S "0.01,0.1,1.0;0.0001,0.001,0.01" --mcmc-samples Outputs/Combined/pack_years.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1

