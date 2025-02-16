library(miamiplot)
library(tidyverse)
library(data.table)
library(CMplot)

#View(gwas_results)
setwd('/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS')
tpy = 'transformed_pack_years_overlap/outfile.fastGWA'
largest_py = 'largest_gwas_pack_years/ErzurumluogluAM_30617275_discovery-stage_meta-analysis_PackYears.txt'
epismoker = 'epismoker_overlap/outfile.fastGWA'
grimage = 'grimage_overlap/outfile.fastGWA'

anno = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno, "V8", "V2")

epismoker_res = read.table(epismoker, header = T)
epismoker_res$STUDY = "Epismoker"
min_smoker_res = epismoker_res[c("SNP", "BETA", "CHR", "POS", "P", "STUDY")]
min_smoker_res = merge(min_smoker_res, anno, by.x="SNP", by.y="V8", all.x = T)
min_smoker_res$SNP = ifelse(is.na(min_smoker_res$V2), min_smoker_res$SNP, min_smoker_res$V2)
min_smoker_res

grimage_res = read.table(grimage, header = T)
grimage_res$STUDY = "grimage"
min_grimage_res = grimage_res[c("SNP", "BETA", "CHR", "POS", "P", "STUDY")]
min_grimage_res = merge(min_grimage_res, anno, by.x="SNP", by.y="V8", all.x = T)

grimage_sig = subset(min_grimage_res, P<5e-8)
write_csv(grimage_sig, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/Grimage_sig_snps.csv")

min_grimage_res$SNP = ifelse(is.na(min_grimage_res$V2), min_grimage_res$SNP, min_grimage_res$V2)
min_grimage_res

tpy = read.table(tpy, header = T)
tpy$STUDY = "Transformed_pack_years"
min_tpy = tpy[c("SNP", "CHR", "BETA", "POS", "P", "STUDY")]
min_tpy = merge(min_tpy, anno, by.x="SNP", by.y="V8", all.x = T)
min_tpy$SNP = ifelse(is.na(min_tpy$V2), min_tpy$SNP, min_tpy$V2)

res = read.table(largest_py, header = T)
res$STUDY = "Pack_years_GWAS"
min_ds = res[c("rsID", "CHROM", "EFFECT_SIZE", "POS", "PVALUE", "STUDY")]
colnames(min_ds) = c("SNP", "CHR", "BETA", "POS", "P", "STUDY")
subset(res, rsID == "rs117836409")

#its all very slow, so I will remove all P vals that are greater than 0.05
min_smoker_res = subset(min_smoker_res, P <= 0.05)
colnames(min_smoker_res) = c("SNP_EpiSmokEr", "BETA_EpiSmokEr", "CHR_EpiSmokEr", "POS_EpiSmokEr", "P_EpiSmokEr", "STUDY_EpiSmokEr", "annot_EpiSmokEr")
head(min_smoker_res)
min_grimage_res = subset(min_grimage_res, P <= 0.05)
colnames(min_grimage_res) = c("SNP_Grimage", "BETA_Grimage", "CHR_Grimage", "POS_Grimage", "P_Grimage", "STUDY_Grimage", "annot_Grimage")
head(min_grimage_res)
min_tpy = subset(min_tpy, P <= 0.05)
colnames(min_tpy) = c("SNP_trans_py", "CHR_trans_py", "BETA_trans_py", "POS_trans_py", "P_trans_py", "STUDY_trans_py", "annot_trans_py")
head(min_tpy)
min_ds = subset(min_ds, P <= 0.05)
colnames(min_ds) = c("SNP_raw", "CHR_raw", "BETA_raw", "POS_raw", "P_raw", "STUDY_raw")
head(min_ds)
subset(min_ds, SNP_raw == "rs117836409")
data = merge(min_smoker_res, min_tpy, by.x="SNP_EpiSmokEr", by.y="SNP_trans_py", all=TRUE)
data = merge(data, min_ds, by.x="SNP_EpiSmokEr", by.y="SNP_raw", all = TRUE)
data = merge(data, min_grimage_res, by.x="SNP_EpiSmokEr", by.y="SNP_Grimage", all = TRUE)

library(dplyr)
data$CHR_raw = as.numeric(data$CHR_raw)
data$Chromosome = coalesce(data$CHR_EpiSmokEr, data$CHR_trans_py, data$CHR_raw, data$CHR_Grimage)
data$Position = coalesce(data$POS_EpiSmokEr, data$POS_trans_py, data$POS_raw, data$POS_Grimage)
dim(data)
df = data.frame("SNP" = data$SNP_EpiSmokEr, 
                "Chromosome" = data$Chromosome,
                "Position" = data$Position,
                "EpiSmoker" = data$P_EpiSmokEr,
                "Transformed_Pack_Years" = data$P_trans_py, 
                "Largest" = data$P_raw,
                "Grimage" = data$P_Grimage
)



saveRDS(df, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/no_filtering.RDS", compress=F)
df = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/four.RDS")
df = subset(df, Chromosome != "X")
SNPs <-  df[
  df$EpiSmoker < 5e-8 |
    df$Transformed_Pack_Years < 5e-8 |
    df$Largest < 5e-8, 1]
SNPs <-  df[
 df$EpiSmoker < 5e-8 |
 df$Transformed_Pack_Years < 5e-8 |
 df$Largest < 5e-8 |
 df$Grimage < 5e-8, 1]

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/plots")

CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=5e-8,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
       file="jpg", file.name="grimage.jpg",
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

colnames(df) = c("CpG", "Chromosome", "Position", "  EpiSmokEr", "  Pack Years", 
                 "Erzurumluoglu et al.", "Raw", "Inverse_transformed" )
CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=c(1e-5, 5e-8),threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col=c("#c94905", "#ffb666"), highlight.col="grey",highlight.cex=0.7,
       file="jpg", file.name="EriSmokEr_Grimage_LargertErz_TransformedPY", ylim=c(0, 18),
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

# Figure 5

setwd('/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS')
df = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/four.RDS")
df = df %>% select(-EpiSmoker)
SNPs <-  df[ 
  !is.na(df$Chromosome) &
  (df$Grimage < 5e-8 |
  df$Transformed_Pack_Years < 5e-8 |
  df$Largest < 5e-8), 1]
SNPs = SNPs[!is.na(SNPs)]
SNPs = SNPs[SNPs != "-"]

colnames(df) = c("CpG", "Chromosome", "Position", "  Pack Years", "Erzurumluoglu et al.", "Grimage Pack Years")
df_05suggestive = df %>% filter(df[,4] < 1e-5 | df[,5] < 1e-5 | df[,6] < 1e-5)
saveRDS(df_05suggestive, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/Figure_5_source_data.RDS")
write.csv(df_05suggestive, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/CMPlot_data/Figure_5_source_data.csv")

CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=c(1e-5, 5e-8),threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col=c("#c94905", "#ffb666"), highlight.col="grey",highlight.cex=0.7,
       file="pdf", file.name="Figure8", ylim=c(0, 18),
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE, lab.cex=1, axis.cex=1.2, legend.cex=1.5)


### Supplemetary Table 10

suggestive <-  subset(df,
 df$Transformed_Pack_Years < 1e-5 |
 df$Largest < 1e-5 |
 df$Grimage < 1e-5)

write.csv(suggestive, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/Grimage_suggestive_snps_full.csv")

results_table = read_csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/Grimage_suggestive_snps_full.csv")
results_table = results_table[-c(1,5)]
df = results_table %>% arrange(Largest, Transformed_Pack_Years, Grimage)

write.csv(df, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/Grimage_suggestive_snps_full.csv")
