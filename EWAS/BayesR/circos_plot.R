library(miamiplot)
library(stringr)
library(dplyr)
library(CMplot)
library(tidyverse)
#View(gwas_results)
path = 'U:/Datastore/IGMM/marioni-lab/Ola/Lab/Smoking/twist_vs_EPIC/Manual_EWAS/'

source = '40_10_EWAS'
twist = readRDS(paste0(path, '40_10_EWAS.RDS'))

twist$study = "TWIST"
twist$chr = str_replace(twist$chr, "chr", "")

# examine non autosomal data
table(twist$chr)
mito_sex = subset(twist, chr %in% c('M', 'X', 'Y'))
table(mito_sex$p < 1e-5) 

# nothing significant, leaving only autosomes
twist$chr = as.numeric(twist$chr)
twist = subset(twist, !is.na(chr))
table(twist$chr)

# reading in EPIC data
epic = read.csv(paste0(path, 'W3_processed.csv'))
epic = epic[-c(1, 10, 11, 12)]
table(epic$chromosome)
epic[epic$chromosome == "chr19_GL949749v2_alt", ] 
# not significant, fine to remove

epic$chromosome = str_replace(epic$chromosome, "chr", "")
epic$chromosome = as.numeric(epic$chromosome)
epic = subset(epic, !is.na(chromosome))

cols = c("cpg", "chr", "pos", "beta", "se", "tstat", "pval", "study")
colnames(twist) = cols
colnames(epic) = cols

# debugging annotation issues
if (FALSE) {
  epic_1 = subset(epic, chr == 1)
  twist_1 = subset(twist, chr == 1)
  one = merge(epic_1, twist_1, by='pos', all = T)
  table(one$chr.x == one$chr.y)
  View(one[1:5000,])
  
  
  link = "../annotations/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed"
  bed <- as.data.frame(read.table(link, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))
  annots_epic = read.table("../annotations/EPIC.hg38.manifest.tsv", header = T)
  
  epic_1 = subset(epic, chr == 1) 
  #watch out! First base in a chromosome is numbered 0 in these coords!!
  
  twist_1 = subset(twist, chr == 1)
  epic_1$pos = epic_1$pos - 1
  one = merge(epic_1, twist_1, by='pos', all = T)
  #table(one$chr.x == one$chr.y) #1492 
  head(one[which(one$chr.x == one$chr.y),])
  
  
  epic_1 = subset(epic, chr == 1)
  twist_1 = subset(twist, chr == 1)
  epic_1$pos = epic_1$pos + 1
  one = merge(epic_1, twist_1, by='pos', all = T)
  #table(one$chr.x == one$chr.y) #44362
  head(one[which(one$chr.x == one$chr.y),])
}


# coordinates in liftover file numbered from 0
epic$pos = epic$pos + 1 
twist$chrpos = paste0(twist$chr, '-', twist$pos)
epic$chrpos = paste0(epic$chr, '-', epic$pos)
all = merge(twist, epic, by="chrpos", all = T)
cols = c("chrpos", "cpg.x", "chr.x", "pos.x", "pval.x", "study.x", 
         "cpg.y", "chr.y", "pos.y", "pval.y", "study.y")
all = all[cols]
table(all$chr.x == all$chr.y)
View(all[1:5000, ])


rm(mito_sex)
rm(epic)
rm(twist)
gc()

# SNP Chromosome Position    trait1     trait2     trait3
# 1 ALGA0000009          1    52297 0.7738187 0.51194318 0.51194318
# 2 ALGA0000014          1    79763 0.7738187 0.51194318 0.51194318
# 3 ALGA0000021          1   209568 0.7583016 0.98405289 0.98405289
# 4 ALGA0000022          1   292758 0.7200305 0.48887140 0.48887140
# 5 ALGA0000046          1   747831 0.9736840 0.22096836 0.22096836
# 6 ALGA0000047          1   761957 0.9174565 0.05753712 0.05753712


gc()
df = data.frame("CpG" = ifelse(is.na(all$cpg.y), all$cpg.x, all$cpg.y), 
                "Chromosome" = ifelse(is.na(all$chr.y), all$chr.x, all$chr.y), 
                "Position" = as.numeric(ifelse(is.na(all$pos.y), all$pos.x, all$pos.y)), 
                "EPIC" = all$pval.y,
                "TWIST" = all$pval.x 
                )


setwd("U:/Datastore/IGMM/marioni-lab/Ola/Lab/Smoking/twist_vs_EPIC/Manual_EWAS/40_10_circos")
saveRDS(df, "40_10_circos_data.RDS", compress=F)
df = subset(df, EPIC < 0.05 | TWIST < 0.05)


CMplot(df, type="p", plot.type="c", r=0.4, col=c("grey30","grey60"),
       chr.labels=paste("Chr",c(1:22), sep=""),
       threshold=c(1e-6,1e-4), cir.chr.h=1.5, amplify=TRUE,
       threshold.lty=c(1,2), threshold.col=c("red", "blue"),
       signal.line=1, signal.col=c("red","green"),
       chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6, outward=FALSE, file="jpg",
       file.name="no_signal_line", dpi=300, file.output=TRUE,
       verbose=TRUE, width=10, height=10)

CMplot(df, type="p", plot.type="c", r=0.4, col=c("grey30","grey60"),
       chr.labels=paste("Chr",c(1:22), sep=""),
       threshold=c(1e-5,1e-3), cir.chr.h=1.5, amplify=TRUE,
       threshold.lty=c(1,2), threshold.col=c("red", "blue"),
       signal.line=1, signal.col=c("red","green"),
       outward=FALSE, file="jpg",
       file.name="no_binning", dpi=300, file.output=TRUE,
       verbose=TRUE, width=10, height=10)

CMplot(df,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.4,
       cir.axis=TRUE, threshold=c(1e-5), threshold.col=c("red"),
       outward=FALSE,cir.axis.col="black",cir.chr.h=1.3,chr.den.col="black",
       file.name="colorful", file="jpg",dpi=300,file.output=TRUE,
       amplify=FALSE, verbose=TRUE,width=10,height=10)

CMplot(df,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.5,
       col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),
       cir.axis=TRUE, threshold=c(1e-5), threshold.col=c("red"),
       outward=FALSE,cir.axis.col="black",cir.chr.h=1.3,chr.den.col="black",
       file.name="threecol", file="jpg",dpi=300,file.output=TRUE,
       amplify=FALSE, verbose=TRUE,width=10,height=10)

# as manhattan

SNPs <-  df[
  df$EPIC < 1e-5 |
  df$TWIST < 1e-5, 1]

saveRDS(df, "40_10_circos_data_renamed.RDS", compress=F)

colnames(df) = c("CpG", "Chromosome", "Position", " EPIC", "TWIST")

CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=c(1e-5, 5e-8),threshold.col=c("grey", "black"),threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col=c("#c94905", "#ffb666"), highlight.col="grey",highlight.cex=0.7,
       file="jpg", file.name="manhattan_two_thresholds", 
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE, trait.legend.pos="left")

CMplot(df,type="p",plot.type="m",LOG10=TRUE,highlight.type="l",highlight=SNPs,
       threshold=c(1e-5),threshold.col=c("black"),threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col=c("#c94905"), highlight.col="grey",highlight.cex=0.7,
       file="jpg", file.name="manhattan", 
       dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE, trait.legend.pos="left")
