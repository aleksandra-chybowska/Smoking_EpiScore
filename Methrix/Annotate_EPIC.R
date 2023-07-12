path = "/Volumes/marioni-lab/Ola/Lab/TWIST/annotations/"
bed_src = "covered_targets_Twist_Methylome_hg38_annotated_separated.bed"
annots_src = "EPIC.hg38.manifest.tsv"
W3_path = "/Volumes/marioni-lab/Ola/Lab/Smoking/twist_vs_EPIC/W3.csv"

bed = read.table(paste0(path, bed_src), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
bed = as.data.frame(bed)
cpgs = read.table(paste0(path, annots_src), header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
cpgs = as.data.frame(cpgs) # 866553
W3 = read.csv(W3_path) # 773860

common_probes = as.data.frame(intersect(W3$cpg, cpgs$Probe_ID)) # 773658
only_in_W3 = setdiff(W3$cpg, cpgs$Probe_ID) # 202  #  finds all rows in x that aren't in y.
only_in_annots = setdiff(cpgs$Probe_ID, W3$cpg) #92895

merged = merge(W3, cpgs, by.x="cpg", by.y="Probe_ID")
merged = merged[c(1:11)]
merged$chromosome = merged$CpG_chrm
merged$position = merged$CpG_beg

write.csv(merged, "/Volumes/marioni-lab/Ola/Lab/Smoking/twist_vs_EPIC/W3_processed.csv")
