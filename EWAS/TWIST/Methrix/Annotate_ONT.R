if (!requireNamespace("annotatr")) {
    BiocManager::install("annotatr", update = FALSE)
}

library("methrix")
library("annotatr")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("tidyverse")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

pass5_40 = readRDS("/Cluster_Filespace/Marioni_Group/Ola/ONT/methrix/5_40_raw_no_snps.RDS")
meth_mat_with_loci <- methrix::get_matrix(m = pass5_40, type = "M", add_loci = TRUE, in_granges = TRUE)

annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic', 'hg38_genes_intronexonboundaries')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)
head(annotations)

dm_annotated = annotate_regions(
  regions = meth_mat_with_loci,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(dm_annotated)

df_dm_annotated = data.frame(dm_annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_dm_annotated))
saveRDS(df_dm_annotated, '/Cluster_Filespace/Marioni_Group/Ola/ONT/methrix/5_40_CpGs_annotated.RDS', compress=F)
ds = readRDS('/Cluster_Filespace/Marioni_Group/Ola/ONT/methrix/5_40_CpGs_annotated.RDS')

annots = df_dm_annotated
annots = ds
ewas = read.csv("/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_full_EWAS_ONT.csv")
head(ewas)

ewas_significant = subset(ewas, ewas$p < 1e-5) # 63
dim(ewas_significant)
annots = annots[c(1:3, 54:63)]
annots = annots %>% unite("cpgs", c("seqnames", "start", "end"), sep="-")
ont_annots = merge(ewas_significant, annots, all.x = TRUE, by.x = "cpg", by.y = "cpgs")
write.csv(ont_annots, "/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_significant_annotated_fully.csv")

# here view produced annots
cleaned_data <- ont_annots %>%
  group_by(cpg) %>%
  summarise(
    chr = first(chr),
    beta = first(beta),
    se = first(se),
    tstat = first(tstat),
    p = first(p),
    annot.symbol = paste(unique(na.omit(annot.symbol)), collapse = ", "),
    annot.type = paste(unique(na.omit(annot.type)), collapse = ", ")
  ) %>% 
  arrange(p)

colnames(cleaned_data) = c("CpG site",	"Chromosome",	"Beta",	"SE", "T-statistic",
  "P",	"Gene symbol",	"Annotation Type"
)
write.csv(cleaned_data, '/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/40_5_46_annotated_hits_cleaned_data.csv', row.names=F)

only_genes_from_annots = annots[!duplicated(annots$annot.symbol), "annot.symbol"]
only_genes_from_annots = na.omit(only_genes_from_annots)

write.table(only_genes_from_annots, '/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_annots_only_genes.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")
saveRDS(annots, '/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_annots.RDS', compress = F)

df = ont_annots
only_genes_from_annots_significant = df[!duplicated(df$annot.symbol), "annot.symbol"]
only_genes_from_annots_significant = na.omit(only_genes_from_annots_significant)

write.table(only_genes_from_annots_significant, '/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_annots_only_genes_significant.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")
write.csv(df, '/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_annotated_hits.csv', row.names=F)

only_genes_from_annots_significant_split = unlist(strsplit(only_genes_from_annots_significant, "-"))
only_genes_from_annots_significant_split = only_genes_from_annots_significant_split[!duplicated(only_genes_from_annots_significant_split)]

write.table(only_genes_from_annots_significant_split, 
  '/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/5_40_46_CpGs_annotated_min_only_genes_significant_split.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")
