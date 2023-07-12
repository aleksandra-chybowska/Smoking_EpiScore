library("methrix")
library("annotatr")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

pass10_40 = readRDS("/Users/shirin/Desktop/10_40_no_snps.RDS")
meth_mat_with_loci <- methrix::get_matrix(m = pass10_40, type = "M", add_loci = TRUE, in_granges = TRUE)

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
saveRDS(df_dm_annotated, '/Volumes/marioni-lab/Ola/Lab/Smoking/TWIST/annotations/40_10_CpGs_annotated.RDS', compress=F)
