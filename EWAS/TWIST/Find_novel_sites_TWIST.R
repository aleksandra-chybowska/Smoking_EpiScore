setwd("C:/Users/s1654019/Documents/Ola/Lab/Smoking/twist_vs_EPIC/Manual_EWAS/")
TWIST_annots = read.csv('40_10_annotated_hits.csv')
df = TWIST_annots[c(1:8, 11, 12, 14)]
df_unique = unique(df)
genes_only = subset(df_unique, !is.na(df_unique$annot.symbol))
novel_1 = subset(genes_only, EWAS_catalogue == "No")

hits_only = readRDS('40_10_EWAS.RDS')
sig_hits = subset(hits_only, p<1e-5)

write.csv(sig_hits, '40_10_not_annotated_hits_wide.csv')

# now try and get full hits with only 
sig_hits$gene = ""
sig_hits$type = ""
sig_hits$in_catalog = ""

for(cpg in sig_hits$cpg) {
  rows = TWIST_annots[which(TWIST_annots$cpg == cpg), ]
  symbols = na.omit(unique(rows$annot.symbol))
  catalog = unique(rows$EWAS_catalogue)
  
  if (length(catalog) != 0) {
    catalog = catalog[catalog != ""]
  }
  
  sig_hits[which(sig_hits$cpg == cpg), "gene"] = paste(symbols, collapse=", ")
  sig_hits[which(sig_hits$cpg == cpg), "type"] = paste(unique(rows$annot.type), collapse=", ")
  sig_hits[which(sig_hits$cpg == cpg), "in_catalog"] = ifelse(length(catalog) == 0, "", catalog)
}
novel = subset(sig_hits, in_catalog == "No")
write.csv(sig_hits, 
          '/Users/shirin/Documents/Edinburgh/Lab/Smoking/twist_vs_EPIC/Manual_EWAS/40_10_annotated_hits_wide.csv', 
          row.names = F)
