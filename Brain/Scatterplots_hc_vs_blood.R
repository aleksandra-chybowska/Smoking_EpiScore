library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(dplyr)

# brain
hc = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/HC_min.RDS")
hc_interesting = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/hc_interesting.R")

hc_results = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/results_hc.RDS")
ahrr = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/AHRR_record_list.R")
ahrr_hc = as.data.frame(ahrr["hc"])
colnames(ahrr_hc) = c("beta", "smoking")

identical(hc_interesting$cg26381592, hc$cg26381592) # TRUE
identical(hc_interesting$cg05575921, ahrr_hc$beta) # TRUE
# OK to use hc_interesting for hc data

# read in blood LBC data (36?)
smoking <- read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
lbc_target <- readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/targets_3489_bloodonly.rds")
rownames(lbc_target) <- lbc_target$Basename
# TODO: Its not mvals but betas, the variable name here is misleading.
lbc_mvals <- readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
lbc_mvals <- t(lbc_mvals)
col <- colnames(hc_interesting)
lbc_mvals_interesting <- as.data.frame(lbc_mvals[,which(colnames(lbc_mvals) %in% col) ])
lbc_mvals_interesting$Sample_Sentrix_ID <- rownames(lbc_mvals_interesting)
head(lbc_mvals_interesting)

ds <- merge(smoking, lbc_target, by.x="lbc36no", by.y="ID")
ds <- merge(ds, lbc_mvals_interesting, by.x="Basename", by.y="Sample_Sentrix_ID")
# ds <- subset(ds, WAVE == 4)
# dim(ds) # 882

# now, subset ds only to interesting sites and keep ID
last_blood = ds[c("ID_raw", "cg05575921", "cg26381592", "WAVE")]
colnames(last_blood) = c("ID_raw", "blood_cg05575921", "blood_cg26381592", "WAVE")

# subset hc to only crucial
only_hc = hc_interesting[c("LBC_ID", "cg05575921", "cg26381592")]
colnames(only_hc) = c("LBC_ID", "hc_cg05575921", "hc_cg26381592")

#LBC360705 - blood sample unavailable

plotting = merge(last_blood, only_hc, by.x = "ID_raw", by.y = "LBC_ID")
plotting <- plotting %>%
  group_by(ID_raw) %>%
  filter(WAVE == max(WAVE)) %>%
  ungroup()

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/")
write.csv(plotting, "Supplemental_Figure_3_Source_Data.csv", row.names = F)

r_cor = cor(x = plotting$blood_cg05575921, y = plotting$hc_cg05575921, method = "pearson")
plot1 = ggplot(plotting, aes(x = blood_cg05575921, y = hc_cg05575921)) +
  geom_point(color = "#0c8bb4") +  # Scatterplot points in blue
  annotate("text", x = 0.4, y = 0.85, 
    label = bquote(r == .(sprintf("%.2f", r_cor))),
    hjust = 0) +
  labs(title = "cg05575921",
       x = NULL,
       y = NULL) +
  theme_minimal()  # Optional theme for a clean look

r_cor = cor(x = plotting$blood_cg26381592, y = plotting$hc_cg26381592, method = "pearson")
plot2 = ggplot(plotting, aes(x = blood_cg26381592, y = hc_cg26381592)) +
  geom_point(color = "#0c8bb4") +  # Scatterplot points in blue
  annotate("text", x = 0, y = 0.04, 
    label = bquote(r == .(sprintf("%.2f", r_cor))),
    hjust = 0) +
  labs(title = "cg26381592",
       x = NULL,
       y = NULL) +
  theme_minimal()  # Optional theme for a clean look
setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/")
pdf("Scatterplot_hc_vs_blood.pdf", width = 8, height = 4)
joined =  plot_grid(plot1, plot2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

y.grob <- textGrob("Brain DNAm", rot=90)

x.grob <- textGrob("Blood DNAm")

#add to plot

grid.arrange(arrangeGrob(joined, left = y.grob, bottom = x.grob))
dev.off()


