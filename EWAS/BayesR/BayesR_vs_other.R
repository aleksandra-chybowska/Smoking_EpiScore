library(tidyverse)

# read in the BayesR results 
setwd("<cluster_home_dir>/TWIST/Annotations")
illumina = read.csv("infinium-methylationepic-v-1-0-b5-manifest-file-from-Illumina.csv")
cat = read.csv('smoking.tsv', sep='\t')
illumina = illumina[c("Name", "CHR_hg38", "Start_hg38", "End_hg38", "Strand_hg38")]
cat_new = merge(cat, illumina, by.y = "Name", by.x = "cpg", all.x = TRUE)
write.csv(cat_new, "smoking_EWAS_catalog_updated.csv")
