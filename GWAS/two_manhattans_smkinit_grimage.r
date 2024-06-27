library(stringr)
library(dplyr)
library(data.table)
library(CMplot)

#View(gwas_results)
setwd('/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS')

# grimage 
# smkinit EUR full
grimage = fread("grimage_overlap/outfile.fastGWA")
smkinit = fread("newest_gwas_published/EUR_without_UKB/EUR_stratified_without_UKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz")
agesmk = fread("newest_gwas_published/EUR_without_UKB/EUR_stratified_without_UKB/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz")
smkces = fread("newest_gwas_published/EUR_without_UKB/EUR_stratified_without_UKB/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz")
cigday = fread("newest_gwas_published/EUR_without_UKB/EUR_stratified_without_UKB/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz")

anno_full = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno_full, "V8", "V2")

grimage$STUDY = "grimage"
min_grimage_res = grimage %>% select("SNP", "BETA", "CHR", "POS", "P", "STUDY")
min_grimage_res = merge(min_grimage_res, anno, by.x="SNP", by.y="V8", all.x = T)
min_grimage_res$SNP = ifelse(is.na(min_grimage_res$V2), min_grimage_res$SNP, min_grimage_res$V2)
min_grimage_res

smkinit$STUDY = "smkinit"
min_smkinit_res = smkinit %>% select("RSID", "BETA", "CHR", "POS", "P", "STUDY")

agesmk$STUDY = "agesmk"
min_agesmk_res = agesmk %>% select("RSID", "BETA", "CHR", "POS", "P", "STUDY")

smkces$STUDY = "smkces"
min_smkces_res = smkces %>% select("RSID", "BETA", "CHR", "POS", "P", "STUDY")

cigday$STUDY = "cigday"
min_cigday_res = cigday %>% select("RSID", "BETA", "CHR", "POS", "P", "STUDY")

min_grimage_res = subset(min_grimage_res, P <= 0.05)
colnames(min_grimage_res) = c("SNP_Grimage", "BETA_Grimage", "CHR_Grimage", "POS_Grimage", "P_Grimage", "STUDY_Grimage", "annot_Grimage")

min_smkinit_res = subset(min_smkinit_res, P <= 0.05)
colnames(min_smkinit_res) = c("SNP_SmkInit", "BETA_SmkInit", "CHR_SmkInit", "POS_SmkInit", "P_SmkInit", "STUDY_SmkInit")

min_agesmk_res = subset(min_agesmk_res, P <= 0.05)
colnames(min_agesmk_res) = c("SNP_agesmk", "BETA_agesmk", "CHR_agesmk", "POS_agesmk", "P_agesmk", "STUDY_agesmk")

min_smkces_res = subset(min_smkces_res, P <= 0.05)
colnames(min_smkces_res) = c("SNP_smkces", "BETA_smkces", "CHR_smkces", "POS_smkces", "P_smkces", "STUDY_smkces")

min_cigday_res = subset(min_cigday_res, P <= 0.05)
colnames(min_cigday_res) = c("SNP_cigday", "BETA_cigday", "CHR_cigday", "POS_cigday", "P_cigday", "STUDY_cigday")

data = merge(min_grimage_res, min_smkinit_res, by.x="SNP_Grimage", by.y="SNP_SmkInit", all=TRUE)
data = merge(data, min_agesmk_res, by.x="SNP_Grimage", by.y="SNP_agesmk", all=TRUE)
data = merge(data, min_smkces_res, by.x="SNP_Grimage", by.y="SNP_smkces", all=TRUE)
data = merge(data, min_cigday_res, by.x="SNP_Grimage", by.y="SNP_cigday", all=TRUE)

head(data)

data$CHR_SmkInit = sub(x = data$CHR_SmkInit, pattern = "chr", replacement = "")
data$CHR_agesmk = sub(x = data$CHR_agesmk, pattern = "chr", replacement = "")
data$CHR_smkces = sub(x = data$CHR_smkces, pattern = "chr", replacement = "")
data$CHR_cigday = sub(x = data$CHR_cigday, pattern = "chr", replacement = "")

data$Chromosome = coalesce(data$CHR_Grimage, 
as.numeric(data$CHR_SmkInit),
as.numeric(data$CHR_agesmk),
as.numeric(data$CHR_smkces),
as.numeric(data$CHR_cigday)
)

data$Position = coalesce(data$POS_Grimage, 
as.numeric(data$POS_SmkInit),
as.numeric(data$POS_agesmk),
as.numeric(data$POS_smkces),
as.numeric(data$POS_cigday)
)

dim(data)

df = data.frame("SNP" = data$SNP_Grimage, 
                "Chromosome" = data$Chromosome,
                "Position" = data$Position,
                "Grimage" = data$P_Grimage,
                "SmkInit" = data$P_SmkInit,
                "AgeSmk" = data$P_agesmk,
                "SmkCes" = data$P_smkces,
                "CigsDay" = data$P_cigday
)

saveRDS(df, "CMPlot_data/newest_plus_grimage_EUR_without_UKB.RDS", compress=F)

df = subset(df, Chromosome != "X")
setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/plots")

SNPs <-  df[
 df$Grimage < 5e-8 |
 df$SmkInit < 5e-8, 1]

CMplot(df,type="p",plot.type="m",LOG10=TRUE,
       threshold=5e-8,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col="red", 
       file="jpg", file.name="grimage_vs_all_EUR_without_UKB",
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

# CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
#        threshold=5e-8,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
#        signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
#        file="jpg", file.name="grimage_vs_all",
#        dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)