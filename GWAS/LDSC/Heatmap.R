library("tidyverse")
library("gplots")
library("R.oo")
library("RColorBrewer")

rg_to_cor_table = function(ds, pheno, gwas) {
  # this function assumes that the ds contains the following columns:
  # pheno_1 - first rg component, eg. EpiSmokEr
  # pheno_2 - second rg component, eg. AgeSmk
  # gwas_type - GWAS type, eg. EUR stratified
  # ds - main data set
  # pheno - "EpiSmokEr" or "Grimage" - first condition
  # gwas - "EUR_stratified"
  
  formatted = ds %>% 
    filter(pheno_1 == pheno & gwas_type == gwas) %>%
    select("pheno_2", "rg") %>%
    rename(!!pheno := rg) %>%
    tibble::column_to_rownames('pheno_2')
  
  return(formatted)
}

bind_rg_tables_by_gwas_type = function(ds, gwas_type) {
  # ds - main data set
  # gwas - GWAS type,eg."EUR_stratified"
  
  unique_phenos = unique(ds$pheno_1)
  
  if (length(unique_phenos) != 2) {
    throw(
      paste0("Only two unique values are supported at this time: ", "pheno_1")
    ) 
  }

  ds1 = rg_to_cor_table(ds, unique_phenos[1], gwas_type)
  ds2 = rg_to_cor_table(ds, unique_phenos[2], gwas_type)
  return(cbind(ds1, ds2))
}

setwd("U:/Ola/Lab/Smoking/LDSC/")

rg_table = read_csv("All_merged_ldsc.csv") %>% 
  rename(pheno_1 = "min_pheno_1", 
         pheno_2 = "pheno 2",
         gwas_type = "min_pheno_2")

glimpse(rg_table)

EUR_full = bind_rg_tables_by_gwas_type(rg_table, "EUR_stratified")
EUR_full = as.matrix(EUR_full)
EUR_without_UKB = bind_rg_tables_by_gwas_type(rg_table, "EUR_stratified_without_UKB")
EUR_without_UKB = as.matrix(EUR_without_UKB)

pdf("plots/EUR_full.pdf", height = 6, width = 6)
heatmap.2(x = EUR_full, 
        dendrogram = "none", 
        cellnote = EUR_full,
        notecol = "darkgray",
        col = brewer.pal(8, "Blues"),
        Rowv = FALSE, 
        Colv = FALSE,
        notecex = 1.5,
        trace = "none",
        key = FALSE,
        cexRow = 1.5,
        cexCol = 1.5,
        margins=c(8,12),
        lwid = c(2,6),
        lhei = c(0.5,8),
        srtCol = 45
        )
dev.off()

pdf("plots/EUR_without_UKB.pdf", height = 6, width = 6)
heatmap.2(x = EUR_without_UKB, 
        dendrogram = "none", 
        cellnote = EUR_without_UKB,
        notecol = "darkgray",
        col = brewer.pal(8, "Blues"),
        Rowv = FALSE, 
        Colv = FALSE,
        notecex = 1.5,
        trace = "none",
        key = FALSE,
        cexRow = 1.5,
        cexCol = 1.5,
        margins=c(8,12),
        lwid = c(2,6),
        lhei = c(0.5,8),
        srtCol = 45
)
dev.off()
