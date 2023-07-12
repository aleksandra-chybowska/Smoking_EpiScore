###########################################################################
##### GENOME-WIDE ASSOCIATION STUDIES ON SMOKING BEHAVIOUR ################
###########################################################################

##################################################################
##### STEP 1. READ IN AND QC PHENOTYPES ##########################
##################################################################

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/")
path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/"
######################################
##### DNA Methylation Data ###########
######################################

#### We want to get IDs of those with DNAm data only and get their FIDs - for file formatting for GCTA ####

# Read in DNAm ids 
samps=readRDS("../GS_methylation/GS20k/GS20k_Targets.rds")
# Read in pedigree data 
ped=read.csv("clinical/pedigree.csv")
# Tidy and combine 
ped=ped[,c(1,2)]
names(ped)=c("FID","IID")
samps1=merge(ped, samps[,c("Sample_Name","Sample_Sentrix_ID")],by.x="IID",by.y="Sample_Name",all.y=T)


#########################################
##### PckYrs + Smoking Status ###########
#########################################

# Read in pack years and smoker status phenotypes
pck=read.csv("updated_smoking_jan_2019/pack_years.csv")
smk=read.csv("updated_smoking_jan_2019/ever_smoke.csv")

## Convert smoking status to three groups - never (1), former (2) and current (3)
smk$new_group=NA
smk[which(smk$ever_smoke %in% 4),"new_group"]=1
smk[which(smk$ever_smoke %in% 2 | smk$ever_smoke %in% 3),"new_group"]=2
smk[which(smk$ever_smoke %in% 1),"new_group"]=3


#########################################
##### Smoking Score #####################
#########################################

# Read in EpiSmokeR score 
epismk=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/epismoker_20k.csv")
# Merge with sample information 
epismk1=merge(samps1,epismk,by.x="Sample_Sentrix_ID",by.y="SampleName")
# Remove NAs 
epismk1=epismk1[which(!is.na(epismk1$smokingScore)),]


#################################
### TOBACCO QUESTIONNAIRES ######
#################################

# Read in tobacco questionnaire information
# tobaccov2=read.csv("PCQ/tobaccov2.csv")
# tobaccov5=read.csv("PCQ/tobaccov2.csv")

#### Here, the tobacco questionnaire was assessed in two batches, therefore we need to harmonise these sub-questionnaires #####
## Take mid-point of ordinal data group in second version (v5) to harmoise to first version (v2)

# tobaccov5$exposure_week2 = NA 
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==0)] = 0
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==1)] = median(1:4)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==2)] = median(5:9)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==3)] = median(10:19)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==4)] = median(20:29)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==5)] = median(30:39)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==6)] = median(40:49)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==7)] = median(50:59)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==8)] = median(60:69)
# tobaccov5$exposure_week2[which(tobaccov5$exposure_week==9)] = 70

# Get environmental tobacco smoking phenotype
# tobaccov5$ets = ifelse(tobaccov5$exposure_week2>=7, 1, 0)
# tobaccov2$ets = ifelse(tobaccov2$exposure_week>=7, 1, 0)



###############################################################
### STEP 2. REMOVE OUTLIERS FOR PACK YEARS + SMOKING SCORE ####
###############################################################

# Define functions to remove outliers 
outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}


# Remove outliers from pack years 
pck$pack_years_clean <- outlierTrim(log(pck$pack_years+1))

# Remove outliers from smoking score 
epismk1$smokingScore <- outlierTrim(epismk1$smokingScore)


###############################################################
##### STEP 3. DEFINE ANALYSIS GROUPS ##########################
###############################################################

## We will compare pack years from (i) former v never smokers 
## and (ii) pack years from never v current smokers 


#### GROUP 1: Former smokers vs never smokers #####

# Extract former + never ids 
form_nev=smk[-which(smk$ever_smoke %in% 1 | smk$ever_smoke %in% 5),]
# Subset pack years data to former and never smokers only 
form_nev_pck=pck[which(pck$Sample_Name %in% form_nev$Sample_Name),]
# Combine with sample information for file formatting 
form_nev_pck=merge(samps1,form_nev_pck,by.x="IID",by.y="Sample_Name")
form_nev_pck=form_nev_pck[!is.na(form_nev_pck$pack_years_clean),]
form_nev_pck1=form_nev_pck[,c("IID","IID")]
names(form_nev_pck1)=c("FID","IID")
form_nev_pck.bayes=form_nev_pck[,c("FID","IID")]
# Save out file 
write.table(form_nev_pck1, paste0(path,  "IDs/former_never_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)
# Save out file for BayesR step - needs FIDs 
write.table(form_nev_pck.bayes, paste0(path, "BayesR/Genetics_Raw/former_never_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)


#### GROUP 2: Current smokers vs never smokers #####

# Extract current + never ids 
cur_nev=smk[which(smk$ever_smoke %in% 1 | smk$ever_smoke %in% 4),]
# > dim(cur_nev)
# [1] 16298     3
# Subset pack years data to current and never smokers only 
cur_nev_pck=pck[which(pck$Sample_Name %in% cur_nev$Sample_Name),]
# Combine with sample information for file formatting 
cur_nev_pck=merge(samps1,cur_nev_pck,by.x="IID",by.y="Sample_Name")
cur_nev_pck=cur_nev_pck[!is.na(cur_nev_pck$pack_years_clean),]
cur_nev_pck1=cur_nev_pck[,c("IID","IID")]
names(cur_nev_pck1)=c("FID","IID")
cur_nev_pck.bayes=cur_nev_pck[,c("FID","IID")]
# Save out file 
write.table(cur_nev_pck1, paste0(path, "IDs/current_never_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)
# Save out file for BayesR step - needs FIDs 
write.table(cur_nev_pck.bayes, paste0(path, "BayesR/Genetics_Raw/current_never_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)



### GROUP 3: Complete pack years information ####

# Remove NAs 
pck1=pck[which(!is.na(pck$pack_years_clean)),]
# Combine with sample information for file formatting 
pck2=merge(samps1,pck1,by.x="IID",by.y="Sample_Name")
pck_ids=pck2[,c("IID","IID")]
names(pck_ids)=c("FID","IID")
pck_ids.bayes=pck2[,c("FID","IID")]
# Save out file 
write.table(pck_ids, paste0(path, "IDs/pack_years_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)
# Save out file for BayesR step - needs FIDs 
write.table(pck_ids.bayes, paste0(path, "BayesR/Genetics_Raw/pack_years_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)



### GROUP 4: Smoking Status ####

# Remove NAs 
smk1=smk[which(!is.na(smk$new_group)),]
# Combine with sample information for file formatting 
smk2=merge(samps1,smk1,by.x="IID",by.y="Sample_Name")
smk_ids=smk2[,c("IID","IID")]
names(smk_ids)=c("FID","IID")
smk_ids.bayes=smk2[,c("FID","IID")]
# Save out file 
write.table(smk_ids, paste0(path, "IDs/smoking_status_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)
# Save out file for BayesR step - needs FIDs 
write.table(smk_ids.bayes, paste0(path, "BayesR/Genetics_Raw/smoking_status_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)


### GROUP 5: EpiSmokeR score ####

# Remove NAs 
epismk2=epismk1[which(!is.na(epismk1$smokingScore)),]
# Combine with sample information for file formatting 
epismk_ids=epismk2[,c("IID","IID")]
names(epismk_ids)=c("FID","IID")
epismk_ids.bayes=epismk2[,c("FID","IID")]
# Save out file 
write.table(epismk_ids, paste0(path, "IDs/smoking_score_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)
# Save out file for BayesR step - needs FIDs 
write.table(epismk_ids.bayes, paste0(path, "BayesR/Genetics_Raw/smoking_score_ids.list"), row.names=F,col.names=F,sep="\t",quote=F)




###############################################################
##### STEP 4. PREAPRE PHENOTYPE FILES #########################
###############################################################

### FILE 1: PACK YEARS ####
pck.phen=pck2[,c("IID","IID","pack_years_clean")]
names(pck.phen)=c("FID", "IID", "phen")
write.table(pck.phen, paste0(path, "Phenotypes/pack_years.phen"), row.names=F, col.names=F, sep =' ', quote = F)

### FILE 2: SMOKING STATUS ####
smk.phen=smk2[,c("IID","IID","new_group")]
names(smk.phen)=c("FID", "IID", "phen")
write.table(smk.phen, paste0(path, "Phenotypes/smoking_status.phen"), row.names=F, col.names=F, sep =' ', quote = F)

### FILE 3: EpiSmokeR SCORE ####
epi.phen=epismk2[,c("IID","IID","smokingScore")]
names(epi.phen)=c("FID", "IID", "phen")
write.table(epi.phen, paste0(path, "Phenotypes/smoking_score.phen"), row.names=F, col.names=F, sep =' ', quote = F)


###############################################################
##### STEP 5. PREAPRE COVARIATE FILES #########################
###############################################################

# Read in covariate file
cov=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/covariates.csv")
# Merge in family IDs for file formatting 
cov=merge(cov,samps1[,c("IID","FID")],by.x="Sample_Name",by.y="IID",all.x=T)
names(cov)[1]="IID"
# Convert sex to number  
cov$sex=ifelse(cov$sex %in% "M", 1, 2)

### FILE 1: BASIC COVARIATES ####
## age and sex 
qcov <- data.frame(FID =cov$IID,
                   IID = cov$IID,
                   age = cov$age)
write.table(qcov, file=paste0(path, "Covariates/quant_basic.qcov"), row.names=F, sep=' ', quote = F)

fcov <- data.frame(FID =cov$IID,
                   IID = cov$IID,
                   sex = cov$sex)
write.table(fcov, file=paste0(path, "Covariates/fact_basic.cov"), row.names=F, sep=' ', quote = F)

### FILE 2: BASIC COVARIATES ####
## age, sex, BMI, genetic PCs  
qcov <- data.frame(FID =cov$IID,
                   IID = cov$IID,
                   age = cov$age,
                   bmi = log(cov$bmi), 
                   PC1 = cov$V3, 
                   PC2 = cov$V4,
                   PC3 = cov$V5,
                   PC4 = cov$V6,
                   PC5 = cov$V7,
                   PC6 = cov$V8,
                   PC7 = cov$V9, 
                   PC8 = cov$V10,
                   PC9 = cov$V11, 
                   PC10 = cov$V12)
write.table(qcov, file=paste0(path, "Covariates/quant_full.qcov"), row.names=F, sep=' ', quote = F)

fcov <- data.frame(FID =cov$IID,
                   IID = cov$IID,
                   sex = cov$sex)
write.table(fcov, file=paste0(path, "Covariates/fact_full.cov"), row.names=F, sep=' ', quote = F)
