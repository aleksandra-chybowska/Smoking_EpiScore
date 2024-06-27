# Draw plots with methylation at similar loci
library(stringr)
library(ggplot2)
library(ggpubr)
home = "/Cluster_Filespace/Marioni_Group/Ola/"

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

## Wave 3
meth_w3 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
meth_w3[which(is.nan(meth_w3))] <- NA
meth_w3[which(is.infinite(meth_w3))] <- NA

# Mean impute missings
w3_nas = apply(meth_w3,1,function(x){which(is.na(x))})
nas = which(unlist(lapply(w3_nas, length))>0)
for(i in names(nas)){
  imp = mean(meth_w3[i,], na.rm=T)
  ind = which(is.na(meth_w3[i,]))
  meth_w3[i,ind] = rep(imp, length(ind))
}

which(is.na(meth_w3))
# numeric(0)
which(is.infinite(meth_w3))
# numeric(0)

# are the sites from TWIST significant?
# cpgs_twist = read.csv(paste0(home, "TWIST/Manual_EWAS/cpgs.csv"))
# this makes no sense since cpgs_twist[i] is a GRange

# list_w3 = read.table("/Cluster_Filespace/Marioni_Group/Daniel/GS_Wave3_EWAS/w3_probes_to_keep.txt")
# meth_w3 = meth_w3[which(row.names(meth_w3) %in% list_w3$V1),] 
meth = t(meth_w3)
meth = m2beta(meth)
rm(meth_w3)
gc()
meth = as.data.frame(meth)

pheno = read.csv(paste0(home, "TWIST/Pheno/2022-07-28-Twist_Sample_Info.csv"))
pheno$smk = ifelse(pheno$ever_smoke == "Never", 0, 1)
table(pheno$smk)
pheno$methrix_id = paste0("X", pheno$Sample_Name)

#data checks
summary(pheno)
# Now subset the DNAm by phenotype file IDs 
meth = meth[which(rownames(meth) %in% pheno$Sample_Sentrix_ID), ]
identical(rownames(meth), pheno$Sample_Sentrix_ID) # FALSE

pheno = pheno[match(rownames(meth), pheno$Sample_Sentrix_ID), ]
identical(rownames(meth), pheno$Sample_Sentrix_ID) # TRUE 


twist = readRDS(paste0(home, "TWIST/Methrix/processed.RDS")) # 40_10
#meth = readRDS(paste0(home, 'TWIST/Methrix/', filename, '.RDS'))
twist = na.omit(twist)
twist = t(twist)
twist = as.data.frame(twist)

# annotation
EPICvsTWIST = read.csv(paste0(home, "/TWIST/Manual_EWAS/EPIC.hg38.manifest.tsv"), sep='\t')

ahrr_twist = twist[c("chr5-373263-373264", 
"chr5-373133-373134", 
"chr5-373135-373136", 
"chr5-373200-373201", 
"chr5-373238-373239",
"chr5-373240-373241",
"chr5-373263-373264"
)]
ahrr_twist$methrix_id = rownames(ahrr_twist)
ahrr_meth = meth["cg05575921"]
ahrr_meth$Sample_Sentrix_ID = rownames(ahrr_meth)

ids = pheno[c("Sample_Sentrix_ID", "methrix_id", "ever_smoke")]

df = merge(ids, ahrr_twist, by.x="methrix_id", by.y = "methrix_id")
df = merge(df, ahrr_meth, by.x="Sample_Sentrix_ID", by.y = "Sample_Sentrix_ID")
cor(df[, "chr5-373263-373264"], df[,"cg05575921"])
cor(df[, "chr5-373133-373134"], df[,"cg05575921"])
cor(df[, "chr5-373135-373136"], df[,"cg05575921"])
cor(df[, "chr5-373200-373201"], df[,"cg05575921"])
cor(df[, "chr5-373238-373239"], df[,"cg05575921"])
cor(df[, "chr5-373240-373241"], df[,"cg05575921"])

df$x = df[,"chr5-373240-373241"]
df$y = df[,"cg05575921"]
pdf("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/plots/AHRR_chr5-373240-373241_ggplot.pdf", 
width = 5, height = 4)
#plot(twist["chr5-373240-373241"], df[,"cg05575921"])
ggplot(df, aes(x=x, y=y, colour=ever_smoke)) +
  geom_point() +
  labs(colour="Smoking status") +
  xlab("TWIST AHRR Methylation (%)") +
  ylab("EPIC AHRR Methylation (%)") +
  geom_smooth(method='lm', colour = "black")
dev.off()

alppl2_twist = twist[c(
"chr2-232419951-232419952",
"chr2-232419961-232419962"
)]

alppl2_twist$methrix_id = rownames(alppl2_twist)
alppl2_meth = meth["cg21566642"]
alppl2_meth$Sample_Sentrix_ID = rownames(alppl2_meth)

ids = pheno[c("Sample_Sentrix_ID", "methrix_id", "ever_smoke")]

df = merge(ids, alppl2_twist, by.x="methrix_id", by.y = "methrix_id")
df = merge(df, alppl2_meth, by.x="Sample_Sentrix_ID", by.y = "Sample_Sentrix_ID")

cor(df[,"chr2-232419951-232419952"], df[,"cg21566642"])
cor(df[,"chr2-232419961-232419962"], df[,"cg21566642"])



pdf("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/plots/ALPPL2_chr2-232419951-232419952.pdf")
plot(df[,"chr2-232419951-232419952"], df[,"cg21566642"])
dev.off()

f2rl3_twist = twist[c(
"chr19-16889774-16889775",
"chr19-16889741-16889742"
)]

f2rl3_twist$methrix_id = rownames(f2rl3_twist)
f2rl3_meth = meth["cg03636183"]
f2rl3_meth$Sample_Sentrix_ID = rownames(f2rl3_meth)

ids = pheno[c("Sample_Sentrix_ID", "methrix_id", "ever_smoke")]

df = merge(ids, f2rl3_twist, by.x="methrix_id", by.y = "methrix_id")
df = merge(df, f2rl3_meth, by.x="Sample_Sentrix_ID", by.y = "Sample_Sentrix_ID")

cor(df[,"chr19-16889774-16889775"], df[,"cg03636183"])
cor(df[,"chr19-16889741-16889742"], df[,"cg03636183"])

pdf("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/plots/F2RL3_chr19-16889774-16889775.pdf")
plot(df[,"chr19-16889774-16889775"], df[,"cg03636183"])
dev.off()

pdf("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/plots/F2RL3_chr19-16889741-16889742.pdf")
plot(df[,"chr19-16889741-16889742"], df[,"cg03636183"])
dev.off()

df$x = df[,"chr19-16889741-16889742"]
df$y = df[,"cg03636183"]
pdf("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/plots/F2RL3_chr19-16889741-16889742-ggplot.pdf", 
width = 5, height = 4)
#plot(twist["chr5-373240-373241"], df[,"cg05575921"])
ggplot(df, aes(x=x, y=y, colour=ever_smoke)) +
  geom_point() +
  labs(colour="Smoking status") +
  xlab("TWIST F2RL3 Methylation (%)") +
  ylab("EPIC F2RL3 Methylation (%)") +
  geom_smooth(method='lm', colour = "black")
dev.off()
df$x = df[,"chr5-373240-373241"]
df$y = df[,"cg05575921"]
pdf("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/plots/AHRR_chr5-373240-373241_ggplot.pdf", 
width = 5, height = 4)
#plot(twist["chr5-373240-373241"], df[,"cg05575921"])
ggplot(df, aes(x=x, y=y, colour=ever_smoke)) +
  geom_point() +
  labs(colour="Smoking status") +
  xlab("TWIST AHRR Methylation (%)") +
  ylab("EPIC AHRR Methylation (%)") +
  geom_smooth(method='lm', colour = "black")
dev.off()