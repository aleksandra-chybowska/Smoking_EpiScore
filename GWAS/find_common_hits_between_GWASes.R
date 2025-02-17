library(tidyverse)
# read in CMPlot prepped file

binary_threshold = function(threshold, col) {
    return(ifelse(col < threshold, 1, 0))
}

setwd("<cluster_home_dir>/Smoking/GWAS/")
gwas = readRDS("CMPlot_data/newest_plus_grimage.RDS")
df = readRDS("<cluster_home_dir>/Smoking/GWAS/CMPlot_data/four.RDS")

t = 5e-8

# Add Erzumogulu here
grimage = subset(gwas, Grimage < t)
largest = subset(df, Largest < t)
gwas_sig = subset(gwas, # nolint
    Grimage < t | SmkInit < t | AgeSmk < t | SmkCes < t | CigsDay < t)



dim(grimage)
write_csv(grimage, "CMPlot_data/hits_overlap_grimage.csv")



gwas_sig[is.na(gwas_sig)] = 1

gwas_sig$Grimage_num = binary_threshold(t, gwas_sig$Grimage)
gwas_sig$SmkInit_num = binary_threshold(t, gwas_sig$SmkInit)
gwas_sig$AgeSmk_num  = binary_threshold(t, gwas_sig$AgeSmk)
gwas_sig$SmkCes_num  = binary_threshold(t, gwas_sig$SmkCes)
gwas_sig$CigsDay_num = binary_threshold(t, gwas_sig$CigsDay)

gwas_sig$shared = gwas_sig$Grimage_num + gwas_sig$SmkInit_num +
    gwas_sig$AgeSmk_num + gwas_sig$SmkCes_num + gwas_sig$CigsDay_num

pdf("plots/newest_plus_grimage_hits_historam.pdf")
hist(gwas_sig$shared)
dev.off()

gwas_sig = subset(gwas_sig, shared > 1)
head(gwas_sig)

pdf("plots/newest_plus_grimage_hits_historam_shared.pdf")
hist(gwas_sig$shared)
dev.off()

# gwas_sig <- gwas_sig %>% select(-contains("_num"))
write.csv(gwas_sig, "CMPlot_data/hits_overlap.csv", row.names = F)

grimage_common = subset(gwas_sig, Grimage_num == 1)
dim(grimage_common)
write_csv(grimage_common, "CMPlot_data/hits_overlap_grimage_common.csv")


### Same for Erzumulogulou

df = readRDS("<cluster_home_dir>/Smoking/GWAS/CMPlot_data/four.RDS")
gwas = df
t = 5e-8

largest_py = subset(df, Largest < t)
also_grimage = subset(largest_py, Grimage < t)
also_grimage$SNP
