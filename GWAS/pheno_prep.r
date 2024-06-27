library(dplyr)

setwd("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/")

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/"
samps = readRDS("../GS_methylation/GS20k/GS20k_Targets.rds")
epismk = read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/epismoker_20k.csv")
grimage = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS_20k_DNAmAge.csv")
grimage = grimage[c("SampleID", "DNAmPACKYRS")]


epismoker = merge(samps, epismk, by.x = "Sample_Sentrix_ID", by.y = "SampleName")
table(is.na(epismoker$smokingScore))
epismoker = subset(epismoker, !is.na(epismoker$smokingScore))

grimage = merge(samps, grimage, by.x = "Sample_Sentrix_ID", by.y = "SampleID")
table(is.na(grimage$DNAmPACKYRS))
grimage = subset(grimage, !is.na(grimage$DNAmPACKYRS))
head(grimage)

py = read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_17865_complete.csv")
min_epismoker = select(epismoker, c("Sample_Sentrix_ID", "smokingScore"))
min_grimage <- select(grimage, c("Sample_Sentrix_ID", "DNAmPACKYRS"))
py$inrt_py = outlierTrim(transform(py$pack_years))

pdf("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/inrt_py.pdf")
hist(py$inrt_py)
dev.off()

pheno <- merge(py, min_epismoker, by.x = "Sample_Sentrix_ID", by.y = "Sample_Sentrix_ID")
pheno <- merge(pheno, min_grimage, by = "Sample_Sentrix_ID")
head(pheno)
dim(pheno)
# 17136     9

PCS = read.csv('/Cluster_Filespace/Marioni_Group/GS/GS_dataset/GS20K_ALL_MAF5_PCA.eigenvec', sep=' ', header = FALSE)
colnames(PCS) = c("FID", "IID", paste0("PC", 1:20))
pheno = merge(PCS, pheno, by.x = "IID", by.y = "Sample_Name",)
# 17130
qcovar = select(pheno, c("IID", "age", 3:22))
qcovar = cbind(qcovar$IID, qcovar)
covar = select(pheno, c("IID", "sex"))
covar = cbind(covar$IID, covar)
covar$sex = as.factor(covar$sex)
pheno_py = select(pheno, c("IID", "pack_years_clean"))
pheno_py = cbind(pheno_py$IID, pheno_py)
pheno_epismoker = select(pheno, c("IID", "smokingScore"))
pheno_epismoker = cbind(pheno_epismoker$IID, pheno_epismoker)
pheno_transformed = select(pheno, c("IID", "inrt_py"))
pheno_transformed = cbind(pheno_transformed$IID, pheno_transformed)
pheno_grimage = select(pheno, c("IID", "DNAmPACKYRS"))
pheno_grimage = cbind(pheno_grimage$IID, pheno_grimage)
raw = select(pheno, c("IID", "pack_years_clean"))
pheno_raw = select(pheno, c("IID", "pack_years"))
pheno_raw = cbind(pheno_raw$IID, pheno_raw)
# check for nas!
table(is.na(qcovar))
table(is.na(covar))
table(is.na(pheno_py))
table(is.na(pheno_epismoker))
table(is.na(pheno_grimage))
table(is.na(pheno_transformed))
pheno_transformed = na.omit(pheno_transformed)
setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/")
write.table(pheno_py, file="transformed_pack_years_overlap/transformed_pack_years_overlap_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(pheno_py[,c("IID", "IID")], file="transformed_pack_years_overlap/samps.txt", sep='\t', col.names=F, row.names=F, quote=F)
write.table(qcovar, file="transformed_pack_years_overlap/qcovar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(covar, file="transformed_pack_years_overlap/covar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)

write.table(pheno_epismoker, file="epismoker_overlap/epismoker_overlap_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(pheno_epismoker[,c("IID", "IID")], file="epismoker_overlap/samps.txt", sep='\t', col.names=F, row.names=F, quote=F)
write.table(qcovar, file="epismoker_overlap/qcovar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(covar, file="epismoker_overlap/covar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)

dir.create("grimage_overlap")
write.table(pheno_grimage, file="grimage_overlap/grimage_overlap_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(pheno_grimage[,c("IID", "IID")], file="grimage_overlap/samps.txt", sep='\t', col.names=F, row.names=F, quote=F)
write.table(qcovar, file="grimage_overlap/qcovar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(covar, file="grimage_overlap/covar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)

dir.create("inrt_overlap")
write.table(pheno_transformed, file="inrt_overlap/inrt_overlap_17124.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(pheno_transformed[,c("IID", "IID")], file="inrt_overlap/samps.txt", sep='\t', col.names=F, row.names=F, quote=F)
write.table(qcovar, file="inrt_overlap/qcovar_17124.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(covar, file="inrt_overlap/covar_17124.phen", sep='\t', row.names=F, quote=F, col.names=F)

dir.create("raw")
write.table(pheno_raw, file="raw/raw_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(pheno_raw[,c("IID", "IID")], file="raw/samps.txt", sep='\t', col.names=F, row.names=F, quote=F)
write.table(qcovar, file="raw/qcovar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(covar, file="raw/covar_17130.phen", sep='\t', row.names=F, quote=F, col.names=F)
