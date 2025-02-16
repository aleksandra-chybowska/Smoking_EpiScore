#read in bed
link = "https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed"
bed = as.data.frame(read.table(link, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed$new_id = paste0("twist", rownames(bed))
saveRDS(bed, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/bed_new_id.RDS')

out <- strsplit(as.character(bed$V4), ',') 
do.call(rbind, out)
require(data.table)

dt = rbindlist(
  lapply(out, function(x) data.table(t(x))),
  fill = TRUE
)

blended = cbind(bed, dt)
saveRDS(blended, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/bed_exploded.RDS')

################ per CpG #####################

library(tidyverse)
library(data.table)

setwd("U:/Ola/Lab/Smoking/TWIST/annotations/")
bed = fread("covered_targets_Twist_Methylome_hg38_annotated_separated.bed")
colnames(bed) = c("chr", "start", "stop", "cpg")
epic = fread("EPIC.hg38.manifest.tsv")
probes = epic$Probe_ID
overlap = intersect(bed$cpg, probes)
saveRDS(overlap, "overlap_between_twist_and_epic.RDS")

hits_only = readRDS('U:/Ola/Lab/Smoking/twist_vs_EPIC/Manual_EWAS/40_10_EWAS.RDS')
sig_hits = subset(hits_only, p<1e-5)

setwd("U:/Ola/Lab/Smoking/twist_vs_EPIC/Manual_EWAS/40_10_circos")
hits = readRDS("40_10_circos_data.RDS")
df = subset(hits, TWIST < 1e-5)

on_illumina = subset(df, CpG %in% overlap) # 10
on_illumina_2 = subset(df, CpG %like% "cg") # 10

common = subset(df, TWIST < 1e-5 & EPIC < 1e-5)
w3 = read.csv("U:/Ola/Lab/Smoking/twist_vs_EPIC/W3_processed.csv")
