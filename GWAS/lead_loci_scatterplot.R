library(tidyverse)
library(ggplot2)

# The overlap of lead loci compared between the GrimAge DNAm pack years GWAS and smoking pack years GWAS 
# can be visualized in scatterplot of the beta coefficients in addition to the current Figure 5

# get GWAS data

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/")
grimage = read.table("grimage_overlap/outfile.fastGWA", header = T)
py = read.table("transformed_pack_years_overlap/outfile.fastGWA", header = T)

# subset(py, P < 5e-8)
# 17:57332520 is 17_57332520_A_G is rs117836409

# find lead loci in FUMA
lead_g = c("2_38244417_G_T", "2_38298139_T_C", "15_78915245_C_T")
lead_py = c("17_57332520_A_G")
lead = c(lead_g, lead_py)

grimage = subset(grimage, SNP %in% lead)
py = subset(py, SNP %in% lead)
gwas_data = data.frame(
    SNP_grimage = grimage$SNP,
    SNP_py = py$SNP,
    Beta_grimage = grimage$BETA,
    Beta_py = py$BETA
)
write.csv(gwas_data, "Supplemental_Figure_5_source_data.csv", row.names=F)
#gwas_data <- read.table("aligned_gwas_betas.txt", header = TRUE)

model <- lm(Beta_grimage ~ Beta_py, data = gwas_data)
r_squared <- summary(model)$r.squared
r_cor <- cor(gwas_data$Beta_grimage, gwas_data$Beta_py,
  method = "pearson", use="complete.obs")


pdf("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/plots/Grimage_DNAm_vs_pack_years_scatterplot_r.pdf",
width = 5, height = 5)
# Plot
# Create the ggplot scatterplot with a line of best fit
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
ggplot(gwas_data, aes(x = Beta_grimage, y = Beta_py)) +
  geom_point(color = "#0c8bb4") +  # Scatterplot points in blue
  annotate("text", x = -1, y = 0.3, 
    label = bquote(r == .(sprintf("%.2f", r_cor))),
    hjust = 0) +
  labs(
    x = "Beta GrimAge DNAm pack years",
    y = "Beta self-reported pack years"
  ) +
  theme_minimal()  # Optional theme for a clean look
dev.off()

# What’s the correlation between smoking pack years
# and GrimAge-based DNAm estimator of smoking pack years in GS
library(data.table)

pack_years = fread("transformed_pack_years_overlap/transformed_pack_years_overlap_17130.phen")
colnames(pack_years) = c("fid", "id", "pack_years")
grimage = fread("grimage_overlap/grimage_overlap_17130.phen")
colnames(grimage) = c("fid", "id", "grimage")
merged = merge(pack_years, grimage, by = "id")
spearman_corr <- cor(merged$pack_years, merged$grimage, method = "spearman")
pearson_corr <- cor(merged$pack_years, merged$grimage, method = "pearson")
print(spearman_corr) #  0.6515853
print(pearson_corr) #  0.703029

