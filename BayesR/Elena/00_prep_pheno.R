library(data.table)
library(dplyr)

setwd("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/")
path <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/"
fam <- "/Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS.fam"
type <- "complete"

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

meanimpute <- function(x) ifelse(is.na(x), mean(x,na.rm=T), x)

ped=read.csv("clinical/pedigree.csv")
ped=ped[,c(1,2)]
names(ped)=c("FID","IID")
target=readRDS("../GS_methylation/GS20k/GS20k_Targets.rds")

pack_years=read.csv("updated_smoking_jan_2019/pack_years.csv")
pack_years=merge(pack_years, target, by.x="Sample_Name", by.y="Sample_Name")

### GROUP 1: Complete pack years information ####
if (type == "complete") {
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years=subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)

}
if (type == "ASSIGN_W1_W3_W4") {
  full = read.csv('/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/data/Full_Dataset.csv') #816 NAs, a bit

  # smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==1 & is.na(smoking_PCQ$cigs_day), smoking_PCQ$packs_day * 20, smoking_PCQ$cigs_day)
  # smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==2, smoking_PCQ$packs_day * 20, smoking_PCQ$cigs_day)
  # smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==3, round_half_up((smoking_PCQ$packs_day * 20) / 2.5), smoking_PCQ$cigarettes_per_day) 
  smoking = data.frame("id" = full$id, "cigs_day" = full$cigs_day, "pack_years" = full$pack_years)
  smoking = merge(smoking, target, by.x = "id", by.y = "Sample_Name")
  smoking = subset(smoking, !is.na(cigs_day))
  sample_size = nrow(smoking)

  pdf(paste0(path, "hist_cigs_day_", sample_size, ".pdf"))
    hist(smoking$cigs_day)
  dev.off()

  smoking$cigs_day_clean = outlierTrim(log(smoking$cigs_day+1))
  smoking=subset(smoking, !is.na(cigs_day_clean))
  sample_size = nrow(packsmoking_years)

  pdf(paste0(path, "hist_log_cigs_day_", sample_size, ".pdf"))
    hist(smoking$cigs_day_clean)
  dev.off()

  out = paste0("cigs_day", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(smoking, paste0(path, out), row.names = F)
}
if (type == "only_W1") {
  pack_years = subset(pack_years, Set=="wave1")
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years=subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}

if (type == "only_W3") {
  pack_years = subset(pack_years, Set=="wave3")
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years=subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}
if (type == "only_W4") {
  pack_years = subset(pack_years, Set=="wave4")
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years=subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}

if (type == "W3_W1") {
  print(type)
  pack_years = subset(pack_years, Set!="wave4")
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years=subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", type, "_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}
if (type == "wave1_wave3_genetic_age_sex") {
  # subset to those with available genetic data
  pack_years = merge(pack_years, ped, by.x="Sample_Name", by.y="IID", all.x = T)
  gen = as.data.frame(fread(fam))
  pack_years = subset(pack_years, Sample_Name %in% gen$V2)

  # subset to wave3, transform pack years
  pack_years = subset(pack_years, Set!="wave4")
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years = subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}

if (type == "wave3_genetic_age_sex") {
  # subset to those with available genetic data
  pack_years = merge(pack_years, ped, by.x="Sample_Name", by.y="IID", all.x = T)
  gen = as.data.frame(fread(fam))
  pack_years = subset(pack_years, Sample_Name %in% gen$V2)

  # subset to wave3, transform pack years
  pack_years = subset(pack_years, Set=="wave3")
  pack_years$pack_years_clean = outlierTrim(log(pack_years$pack_years+1))
  pack_years = subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "hist_pack_years_clean_", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}

if (type == "no_transform") {
  pack_years$pack_years = outlierTrim(pack_years$pack_years)
  pack_years=subset(pack_years, !is.na(pack_years))
  sample_size = nrow(pack_years)

  pdf(paste0(path, "plots/hist_pack_years_no_transform", nrow(pack_years), ".pdf"))
    hist(pack_years$pack_years)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)

}

if (type == "complete_genetic") {
  pack_years=merge(pack_years, ped, by.x="Sample_Name", by.y="IID", all.x = T)
  gen=as.data.frame(fread(fam))
  pack_years = subset(pack_years, Sample_Name %in% gen$V2)

  pack_years$pack_years_clean <- outlierTrim(log(pack_years$pack_years+1))
  pack_years=subset(pack_years, !is.na(pack_years_clean))
  sample_size = nrow(pack_years)

  # subset to those with genetic data
  pdf(paste0(path, "hist_pack_years_", sample_size, ".pdf"))
    hist(pack_years$pack_years)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(pack_years, paste0(path, out), row.names = F)
}

if (type == "former_current") {
  smk <- read.csv("updated_smoking_jan_2019/ever_smoke.csv")

  ## Convert smoking status to three groups - never (1), former (2) and current (3)
  smk$new_group=NA
  smk[which(smk$ever_smoke %in% 4),"new_group"]=1
  smk[which(smk$ever_smoke %in% 2 | smk$ever_smoke %in% 3),"new_group"]=2
  smk[which(smk$ever_smoke %in% 1),"new_group"]=3
  former_current = subset(smk, new_group %in% c(2, 3))

  ds <- merge(pack_years, former_current, by="Sample_Name")
  ds$pack_years_clean <- outlierTrim(log(ds$pack_years+1))
  ds <- subset(ds, !is.na(pack_years_clean))
  sample_size = nrow(ds)

  pdf(paste0(path, "plots/former_current_pack_years.pdf"))
  hist(ds$pack_years_clean)
  dev.off()

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(ds, paste0(path, out), row.names = F)
}

if (type == "bmi") {
  ### BMI for sensitivity analysis
  cov=read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/common_covariates.csv")
  cov=cov[c(1,3)]

  ds=merge(pack_years, cov, by="Sample_Name")
  summary(ds$bmi) #109 NAs
  meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  ds$bmi=meanimpute(ds$bmi)
  summary(ds$bmi) #0 NAs

  pdf(paste0(path, "hist_bmi.pdf"))
  hist(log(ds$bmi))
  dev.off()

  ds$bmi=log(ds$bmi)
  sample_size = nrow(ds)

  out = paste0("pack_years_", sample_size, "_", type, ".csv")
  print(paste("Prepped pheno file: ", out))
  write.csv(ds, paste0(path, out), row.names = F)
}