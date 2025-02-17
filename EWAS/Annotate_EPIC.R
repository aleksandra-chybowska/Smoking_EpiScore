path = "<cluster_home_dir>/TWIST/Annotations/"
bed_src = "covered_targets_Twist_Methylome_hg38_annotated_separated.bed"
annots_src = "infinium-methylationepic-v-1-0-b5-manifest-file-from-Illumina.csv"
W3_path = "<cluster_home_dir>/TWIST/Manual_EWAS/W3_46.csv"

# path = "<datastore_home_dir>/Lab/Smoking/TWIST/annotations/"
# bed_src = "covered_targets_Twist_Methylome_hg38_annotated_separated.bed"
# annots_src = "EPIC.hg38.manifest.tsv"
# W3_path = "<datastore_home_dir>/Lab/Smoking/twist_vs_EPIC/W3.csv"

bed = read.table(paste0(path, bed_src), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
bed = as.data.frame(bed)
cpgs = read.csv(paste0(path, annots_src))
W3 = read.csv(W3_path) # 773860

common_probes = as.data.frame(intersect(W3$cpg, cpgs$Name)) # 773658 # 773658
only_in_W3 = setdiff(W3$cpg, cpgs$Name) # 202  #  finds all rows in x that aren't in y.
only_in_annots = setdiff(cpgs$Name, W3$cpg) #92895 # 92276

merged = merge(W3, cpgs, by.x="cpg", by.y="Name")
merged = cbind(merged[1:8], merged[56:58])
merged$chromosome = merged$CHR_hg38
merged$position = merged$Start_hg38

write.csv(merged, "<cluster_home_dir>/TWIST/Manual_EWAS/W3_46_processed_new_annots.csv")

### for FUMA

extract_gene_id = function(all_genes) {
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

  all_genes <- sapply(all_genes, function(x) {
    # Split each element by ";", remove duplicates with unique()
    unique(strsplit(x, ";")[[1]])
  })
  head(all_genes)
  all_genes <- unlist(unique(all_genes))
  return(all_genes)
}

setwd("<cluster_home_dir>/TWIST/Annotations/")
otherEPIC <- read.csv("infinium-methylationepic-v-1-0-b5-manifest-file-from-Illumina.csv")
anno <- otherEPIC[c("CHR_hg38","Name", "Start_hg38","UCSC_RefGene_Name", "Strand_hg38")]

w3 <- read.csv("<cluster_home_dir>/TWIST/Manual_EWAS/W3.csv") 
head(w3)

ewas_significant = subset(ewas, ewas$pval < 1e-5) # 97
epic = merge(ewas_significant, anno, by.x = "cpg", by.y = "Name")

background = extract_gene_id(otherEPIC$UCSC_RefGene_Name)
front = extract_gene_id(epic$UCSC_RefGene_Name)

write.table(front, 
  '../Manual_EWAS/epic/EPIC_FUMA_front.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")

write.table(background, 
  '../Manual_EWAS/epic/EPIC_FUMA_background.csv', 
  row.names = F, col.names = F, quote = FALSE, sep=",")
