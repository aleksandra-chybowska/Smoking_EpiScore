# Load necessary library
library(dplyr)
library(readxl)


setwd("<local_home_dir>/Desktop/")
# Read in the data (assuming your file is named 'data.csv' and is comma-separated)

cleaned_data_ont = read.csv("40_5_46_annotated_hits_cleaned_data.csv")
# Display the first few rows of the cleaned data
ont = cleaned_data_ont
twist =  read.csv("40_10_46_annotated_hits_cleaned_data.csv")
head(twist)
twist_genes = twist$Gene.symbol
ont_genes = ont$Gene.symbol

twist_genes <- twist_genes[twist_genes != ""]
twist_genes <- unlist(strsplit(twist_genes, ", "))
ont_genes <- ont_genes[ont_genes != ""]
ont_genes <- unlist(strsplit(ont_genes, ", "))


# Display the cleaned list of genes

cat = read.csv('smoking_EWAS_catalog_updated.csv') # as I need coords in hg38

twist_genes_in_cat_name = twist_genes[twist_genes %in% cat$gene]
twist_genes_not_in_cat_name = twist_genes[!(twist_genes %in% cat$gene)]
# first, do it by gene 

twist <- transform(twist, 
                chr = sub("-.*", "", CpG.site),
                start = as.numeric(sub(".*-(\\d+)-\\d+$", "\\1", CpG.site)) - 1, # gotta do it because of the way epic is indexed
                stop = as.numeric(sub(".*-(\\d+)$", "\\1", CpG.site))
)
twist_sig = subset(twist, P<3.6e-8)
twist_sug = subset(twist, P<1e-5)

twist_loci_in_cat = subset(twist_sig, chr %in% cat$CHR_hg38 & start %in% cat$Start_hg38) # 2
twist_loci_not_in_cat = subset(twist_sig, !(chr %in% cat$CHR_hg38 & start %in% cat$Start_hg38)) # 31

twist_sug_loci_in_cat = subset(twist_sug, chr %in% cat$CHR_hg38 & start %in% cat$Start_hg38) # 4
twist_sug_loci_not_in_cat = subset(twist_sug, !(chr %in% cat$CHR_hg38 & start %in% cat$Start_hg38)) # 98

# overlap twist epic

test = c("cg14391737", "cg03636183", "cg21566642", "cg05575921")
twist_by_colname = subset(twist, start %in% 373263)
cat_by_colname = subset(cat, pos %in% 373265)


ont <- transform(ont, 
                 chr1 = sub("-.*", "", CpG.site),
                 start1 = as.numeric(sub(".*-(\\d+)-\\d+$", "\\1", CpG.site)) - 1,
                 stop1 = as.numeric(sub(".*-(\\d+)$", "\\1", CpG.site))
)

# Genome wide significant 
ont_sig = subset(ont, P<3.6e-8)
ont_loci_in_cat = subset(ont_sig, chr1 %in% cat$CHR_hg38 & start1 %in% cat$Start_hg38) # 1
ont_loci_not_in_cat = subset(ont_sig, !(chr1 %in% cat$CHR_hg38 & start1 %in% cat$Start_hg38)) # 62
# suggestive
ont_sug = subset(ont, P<1e-5)
ont_loci_in_cat = subset(ont_sug, chr1 %in% cat$CHR_hg38 & start1 %in% cat$Start_hg38) # 1
ont_loci_not_in_cat = subset(ont_sug, !(chr1 %in% cat$CHR_hg38 & start1 %in% cat$Start_hg38)) # 62

ont_by_colname = subset(ont, start1 %in% 373262)

bayesR = read.csv("BayesR.csv")
bayesR_in_cat = subset(bayesR, !is.na(EWAS.catalog.n))
bayesR_in_cat_high_P = subset(bayesR_in_cat, EWAS.catalog.P < 1e-7)

novel_80 = subset(bayesR, is.na(EWAS.catalog.n))
