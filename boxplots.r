library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(gridExtra)

hc = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/HC_min.RDS")
hc_results = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/results_hc.RDS")
ahrr = readRDS("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/AHRR_record_list.R")

ahrr_hc = as.data.frame(ahrr["hc"])
colnames(ahrr_hc) = c("beta", "smoking")
# Brain 
hc$smoking_cat = factor(hc$smoking_cat, labels=c("Never", "Ex", "Current"))
ahrr_hc$smoking = factor(ahrr_hc$smoking, labels=c("Never", "Ex", "Current"))
pdf(file="/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/Plots/AHRR_hc_boxplot_2.pdf")
par(mfrow=c(1,2))
p1 = ggplot(hc, aes(x=as.factor(smoking_cat), y=cg26381592, fill=smoking_cat)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values=c("#102033", "#285686", "#47A0F4")) +
  geom_jitter(color="black", size=0.4, alpha=0.9, show.legend = FALSE,) +
  xlab("") +
  ylab("")
  # cg26381592

p2 = ggplot(ahrr_hc, aes(x=as.factor(smoking), y=beta, fill=smoking)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values=c("#102033", "#285686", "#47A0F4")) +
  geom_jitter(color="black", size=0.4, alpha=0.9, show.legend = FALSE) +
  xlab("") +
  ylab("")
  # AHRR
dev.off()



col <- c("cg05575921", colnames(hc))
# Blood

smoking <- read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
lbc_target <- readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/targets_3489_bloodonly.rds")
rownames(lbc_target) <- lbc_target$Basename
lbc_mvals <- readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
lbc_mvals <- t(lbc_mvals)
lbc_mvals_interesting <- as.data.frame(lbc_mvals[,which(colnames(lbc_mvals) %in% col) ])
lbc_mvals_interesting$Sample_Sentrix_ID <- rownames(lbc_mvals_interesting)

# rm(lbc_mvals)
# gc()

ds <- merge(smoking, lbc_target, by.x="lbc36no", by.y="ID")
ds <- merge(ds, lbc_mvals_interesting, by.x="Basename", by.y="Sample_Sentrix_ID")

head(ds[c("cg05575921", "cg26381592")])

pdf(file="/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/Plots/AHRR_blood_boxplot.pdf")
par(mfrow=c(1,2))
p3 = ggplot(ds, aes(x=as.factor(smokcat_w1), y=cg05575921, fill=as.factor(smokcat_w1))) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values=c("#102033", "#285686", "#47A0F4")) +
  geom_jitter(color="black", size=0.4, alpha=0.9, show.legend = FALSE) +
  xlab("") +
  ylab("AHRR")
 # ylab("")

p4 = ggplot(ds, aes(x=as.factor(smokcat_w1), y=cg26381592, fill=as.factor(smokcat_w1))) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values=c("#102033", "#285686", "#47A0F4")) +
  geom_jitter(color="black", size=0.4, alpha=0.9, show.legend = FALSE) +
  xlab("") +
  ylab("cg26381592")
  # 
dev.off()

pdf(file="/Cluster_Filespace/Marioni_Group/Ola/Smoking/Brain_vs_Blood/Plots/Grid_jitter.pdf", width = 4, height = 4)
grid.arrange(p3, p2, p4, p1,  ncol=2, nrow = 2)
dev.off()