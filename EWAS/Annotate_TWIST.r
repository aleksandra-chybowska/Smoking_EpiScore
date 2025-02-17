library(tidyverse)

annots = readRDS('<cluster_home_dir>/TWIST/Annotations/40_10_CpGs_annotated.RDS')
ewas = read.csv('<cluster_home_dir>/TWIST/Manual_EWAS/pass10_40_46_full_EWAS_TWIST.csv')
head(ewas)


annots = annots[c(1:3, 54:63)]
annots = annots %>% unite("cpgs", c("seqnames", "start", "end"), sep="-")
saveRDS(annots, '<cluster_home_dir>/TWIST/Annotations/40_10_CpGs_annotated_min.RDS', compress = F)


ewas_significant = subset(ewas, ewas$p < 1e-5) # 102
df = merge(ewas_significant, annots, by.x = "cpg", by.y="cpgs")
write.csv(df, '<cluster_home_dir>/TWIST/Annotations/40_10_46_annotated_hits.csv', row.names=F)

cleaned_data <- df %>%
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
#prepare supplemental table 5
write.csv(cleaned_data, '<cluster_home_dir>/TWIST/Annotations/40_10_46_annotated_hits_cleaned_data.csv', row.names=F)


only_genes_from_annots = annots[!duplicated(annots$annot.symbol), "annot.symbol"]
only_genes_from_annots = na.omit(only_genes_from_annots)
write.table(only_genes_from_annots, 
  '<cluster_home_dir>/TWIST/Annotations/40_10_46samples_CpGs_annotated_min_only_genes.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")

only_genes_from_annots_significant = df[!duplicated(df$annot.symbol), "annot.symbol"]
only_genes_from_annots_significant = na.omit(only_genes_from_annots_significant)

write.table(only_genes_from_annots_significant, 
  '<cluster_home_dir>/TWIST/Annotations/40_10_CpGs_annotated_min_only_genes_significant.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")

only_genes_from_annots_significant_split = unlist(strsplit(only_genes_from_annots_significant, "-"))
only_genes_from_annots_significant_split = only_genes_from_annots_significant_split[!duplicated(only_genes_from_annots_significant_split)]

write.table(only_genes_from_annots_significant_split, 
  '<cluster_home_dir>/TWIST/Annotations/40_10_CpGs_annotated_min_only_genes_significant_split.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")

############
# Now EPIC

anno <- readRDS("<filespace_marioni_group_dir>/Daniel/EPIC_AnnotationObject_df.rds")
anno <- anno[c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
ewas <- read.csv('<cluster_home_dir>/TWIST/Manual_EWAS/W3_46.csv')

head(ewas)
ewas_significant = subset(ewas, ewas$pval < 1e-5) #44
epic = merge(ewas_significant, anno, by.x = "cpg", by.y = "Name")
write.csv(epic, '<cluster_home_dir>/TWIST/Manual_EWAS/epic/40_10_46_annotated_hits.csv', row.names=F)
