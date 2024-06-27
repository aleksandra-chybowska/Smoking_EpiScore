library(miamiplot)
library(stringr)
library(dplyr)
library(data.table)
library(CMplot)

#View(gwas_results)
setwd('/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS')
tpy = 'transformed_pack_years_overlap/outfile.fastGWA'
largest_py = 'largest_gwas_pack_years/ErzurumluogluAM_30617275_discovery-stage_meta-analysis_PackYears.txt'
epismoker = 'epismoker_overlap/outfile.fastGWA'
raw_py = 'raw/outfile.fastGWA'
inrt_py = 'inrt_overlap/outfile.fastGWA'

anno_full = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno, "V8", "V2")

epismoker_res = read.table(epismoker, header = T)
epismoker_res$STUDY = "Epismoker"
min_smoker_res = epismoker_res[c("SNP", "BETA", "CHR", "POS", "P", "STUDY")]
min_smoker_res = merge(min_smoker_res, anno, by.x="SNP", by.y="V8", all.x = T)
min_smoker_res$SNP = ifelse(is.na(min_smoker_res$V2), min_smoker_res$SNP, min_smoker_res$V2)
min_smoker_res

tpy = read.table(tpy, header = T)
tpy$STUDY = "Transformed_pack_years"
min_tpy = tpy[c("SNP", "CHR", "BETA", "POS", "P", "STUDY")]
min_tpy = merge(min_tpy, anno, by.x="SNP", by.y="V8", all.x = T)
min_tpy$SNP = ifelse(is.na(min_tpy$V2), min_tpy$SNP, min_tpy$V2)

raw_py = read.table(raw_py, header = T)
raw_py$STUDY = "Raw_pack_years"
min_raw = raw_py[c("SNP", "CHR", "BETA", "POS", "P", "STUDY")]
min_raw = merge(min_raw, anno, by.x="SNP", by.y="V8", all.x = T)
min_raw$SNP = ifelse(is.na(min_raw$V2), min_raw$SNP, min_raw$V2)

inrt_py = read.table(inrt_py, header = T)
inrt_py$ = "Inverse_transformed_pack_years"
min_inrt = inrt_py[c("SNP", "CHR", "BETA", "POS", "P", "STUDY")]
min_inrt = merge(min_inrt, anno, by.x="SNP", by.y="V8", all.x = T)
min_inrt$SNP = ifelse(is.na(min_inrt$V2), min_inrt$SNP, min_inrt$V2)

res = read.table(largest_py, header = T)
res$STUDY = "Pack_years_GWAS"
min_ds = res[c("rsID", "CHROM", "EFFECT_SIZE", "POS", "PVALUE", "STUDY")]
colnames(min_ds) = c("SNP", "CHR", "BETA", "POS", "P", "STUDY")

#its all very slow, so I will remove all P vals that are greater than 0.05
min_smoker_res = subset(min_smoker_res, P <= 0.05)
head(min_smoker_res)
colnames(min_smoker_res) = c("SNP_EpiSmokEr", "BETA_EpiSmokEr", "CHR_EpiSmokEr", "POS_EpiSmokEr", "P_EpiSmokEr", "STUDY_EpiSmokEr", "annot_EpiSmokEr")
min_tpy = subset(min_tpy, P <= 0.05)
head(min_tpy)
colnames(min_tpy) = c("SNP_trans_py", "CHR_trans_py", "BETA_trans_py", "POS_trans_py", "P_trans_py", "STUDY_trans_py", "annot_trans_py")
min_ds = subset(min_ds, P <= 0.05)
head(min_ds)
min_raw = subset(min_raw, P <= 0.05)
head(min_raw)
colnames(min_raw) = c("SNP_raw", "CHR_raw", "BETA_raw", "POS_raw", "P_raw", "STUDY_raw", "annot_raw")
min_inrt = subset(min_inrt, P <= 0.05)
head(min_inrt)
colnames(min_inrt) = c("SNP_inrt", "CHR_inrt", "BETA_inrt", "POS_inrt", "P_inrt", "STUDY_inrt", "annot_inrt")

data = merge(min_smoker_res, min_tpy, by.x="SNP_EpiSmokEr", by.y="SNP_trans_py", all=TRUE)
data = merge(data, min_ds, by.x="SNP_EpiSmokEr", by.y="SNP", all = TRUE)
data = merge(data, min_raw, by.x="SNP_EpiSmokEr", by.y="SNP_raw", all = TRUE)
data = merge(data, min_inrt, by.x="SNP_EpiSmokEr", by.y="SNP_inrt", all = TRUE)
head(data)

library(dplyr)
data$CHR = as.numeric(data$CHR)
data$Chromosome = coalesce(data$CHR_EpiSmokEr, data$CHR_trans_py, data$CHR, data$CHR_raw, data$CHR_inrt)
data$Position = coalesce(data$POS_EpiSmokEr, data$POS_trans_py, data$POS, data$POS_raw, data$POS_inrt)
dim(data)
df = data.frame("SNP" = data$SNP_EpiSmokEr, 
                "Chromosome" = data$Chromosome,
                "Position" = data$Position,
                "EpiSmoker" = data$P_EpiSmokEr,
                "Transformed_Pack_Years" = data$P_trans_py, 
                "Largest" = data$P,
                "Raw" = data$P_raw,
                "Inverse_transformed" = data$P_inrt
)

saveRDS(df, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/five.RDS", compress=F)
df = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/five.RDS")
df = subset(df, Chromosome != "X")
setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/plots")

df2 = df[c("SNP", "Chromosome", "Position", "EpiSmoker", "Largest", "Inverse_transformed")]
SNPs <-  df2[
 df$EpiSmoker < 5e-8 |
 df$Largest < 5e-8 |
 df$Inverse_transformed < 5e-8 , 1]

CMplot(df2,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=5e-8,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
       file="jpg", file.name="epi_largest_inverse_transformed.jpg",
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

table(df$EpiSmoker < 1e-15)
df[which(df$EpiSmoker < 1e-15), "EpiSmoker"] = 1e-15
table(df$Transformed_Pack_Years < 1e-15)
table(df$Largest < 1e-15)
df[which(df$Largest < 1e-15), "Largest"] = 1e-15

SNPs <-  df[
  df$EpiSmoker < 5e-8 |
  df$Transformed_Pack_Years < 5e-8 |
  df$Largest < 5e-8, 1]

colnames(df) = c("CpG", "Chromosome", "Position", "  EpiSmokEr", "  Pack Years", "Erzurumluoglu et al." )
CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=c(1e-5, 5e-8),threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col=c("#c94905", "#ffb666"), highlight.col="grey",highlight.cex=0.7,
       file="jpg", file.name="final", ylim=c(0, 18),
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

# SNP Chromosome Position    trait1     trait2     trait3
# 1 ALGA0000009          1    52297 0.7738187 0.51194318 0.51194318
# 2 ALGA0000014          1    79763 0.7738187 0.51194318 0.51194318
# 3 ALGA0000021          1   209568 0.7583016 0.98405289 0.98405289
# 4 ALGA0000022          1   292758 0.7200305 0.48887140 0.48887140
# 5 ALGA0000046          1   747831 0.9736840 0.22096836 0.22096836
# 6 ALGA0000047          1   761957 0.9174565 0.05753712 0.05753712

significant = subset(df, EpiSmoker < 5e-8 | Largest < 5e-8 | Inverse_transformed < 5e-8 | Raw < 5e-8 | Transformed_Pack_Years < 5e-8)
write.csv(significant, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/three/five_significant.csv")
