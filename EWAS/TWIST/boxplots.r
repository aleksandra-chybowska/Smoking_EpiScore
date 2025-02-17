library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(gridExtra)

setwd("<cluster_home_dir>/Smoking/Brain_vs_Blood")
hc = readRDS("<cluster_home_dir>/Smoking/Brain_vs_Blood/HC_min.RDS")
# hippocampal methylation values
hc_results = readRDS("<cluster_home_dir>/Smoking/Brain_vs_Blood/results_hc.RDS")
ahrr = readRDS("<cluster_home_dir>/Smoking/Brain_vs_Blood/AHRR_record_list.R")

ahrr_hc = as.data.frame(ahrr["hc"])
colnames(ahrr_hc) = c("beta", "smoking")
# Brain 
hc$smoking_cat = factor(hc$smoking_cat, labels=c("Never", "Ex", "Current"))
ahrr_hc$smoking = factor(ahrr_hc$smoking, labels=c("Never", "Ex", "Current"))

p1 = ggplot(hc, aes(x=as.factor(smoking_cat), y=cg26381592, color=smoking_cat)) + # fill=name allow to automatically dedicate a color for each group
  scale_color_manual(values=c("#C94905", "#218F60", "#155DA2")) +
  geom_point(size=0.7, alpha=0.9,show.legend = FALSE, position=position_jitter(w = 0.15, h = 0)) +
  xlab("") +
  ylab("") +
  ylim(0, 0.04)
  # cg26381592

p2 = ggplot(ahrr_hc, aes(x=as.factor(smoking), y=beta, color=smoking)) + # fill=name allow to automatically dedicate a color for each group
  scale_color_manual(values=c("#C94905", "#218F60", "#155DA2")) +
  geom_point(size=0.7, alpha=0.9,show.legend = FALSE, position=position_jitter(w = 0.15, h = 0)) +
  xlab("") +
  ylab("") +
  ylim(0.4, 1)
  # AHRR
brain_other = hc[c("cg26381592", "smoking_cat")]

write.csv(brain_other, "Figure4_source_data_other_HC.csv", row.names = F)
write.csv(ahrr_hc, "Figure4_source_data_AHRR_HC.csv", row.names = F)
col <- c("cg05575921", colnames(hc))

smoking <- read.csv("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/LBC/pheno_min_1072.csv")
lbc_target <- readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/targets_3489_bloodonly.rds")
rownames(lbc_target) <- lbc_target$Basename
# TODO: Its not mvals but betas, the variable name here is misleading.
lbc_mvals <- readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
lbc_mvals <- t(lbc_mvals)
lbc_mvals_interesting <- as.data.frame(lbc_mvals[,which(colnames(lbc_mvals) %in% col) ])
lbc_mvals_interesting$Sample_Sentrix_ID <- rownames(lbc_mvals_interesting)

ds <- merge(smoking, lbc_target, by.x="lbc36no", by.y="ID")
ds <- merge(ds, lbc_mvals_interesting, by.x="Basename", by.y="Sample_Sentrix_ID")
ds <- subset(ds, WAVE == 1)
dim(ds) # 882

ds$smoking = factor(ds$smokcat_w1, labels=c("Never", "Ex", "Current"))
head(ds[c("cg05575921", "cg26381592", "smoking")])
blood = ds[c("cg05575921", "cg26381592", "smoking")]
write.csv(blood, "Figure4_source_data_Blood.csv", row.names = F)

p3 = ggplot(ds, aes(x=as.factor(smoking), y=cg05575921, fill=as.factor(smoking))) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values=c("#C94905", "#218F60", "#155DA2")) +
  # geom_jitter(color="black", size=0.4, alpha=0.9, show.legend = FALSE) +
  xlab("") +
  ylab("cg05575921")

p4 = ggplot(ds, aes(x=as.factor(smoking), y=cg26381592, fill=as.factor(smoking))) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_manual(values=c("#C94905", "#218F60", "#155DA2")) +
  xlab("") +
  ylab("cg26381592") +
  ylim(0, 0.04)

pdf(file="<cluster_home_dir>/Smoking/Brain_vs_Blood/Plots/Grid.pdf", width = 4, height = 4)

grid.arrange(arrangeGrob(p3, p4, top = grid::textGrob("Blood", x = 0.5, hjust = 0)),
             arrangeGrob(p2, p1, top = grid::textGrob("Hippocampus", x = 0.35, hjust = 0)), ncol = 2)
dev.off()