library(biomaRt)
library(tidyverse)

listMarts()
ensembl <- useMart("ENSEMBL_MART_SNP")
searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl <- useDataset(dataset = "hsapiens_snp", mart = ensembl)
genes <- useMart("ENSEMBL_MART_ENSEMBL")
genes <- useDataset(dataset = "hsapiens_gene_ensembl", mart = genes)
  
setwd("<datastore_home_dir>/Lab/Smoking/LDSC/overlap_between_hits/")
gwas <- read_csv("hits_overlap_grimage.csv")



rsids <- gwas$SNP
# overlap with erzumoglou
rsids = c("rs569207", "rs578776", "rs6495308", "rs6495309", "rs667282", "rs8042374", "rs938682" )
# overlap with cigs day newest
rsids = c("rs11072774","rs11637630","rs117832812","rs12440014","rs12593950","rs12910984","rs13329271","rs1700006","rs1711731","rs2456020","rs28395178","rs28534575","rs28669908","rs28681284","rs3743078","rs376499689","rs3825845","rs4887069","rs4887072","rs518425","rs569207","rs576982","rs578776","rs637137","rs6495308","rs6495309","rs667282","rs684513","rs7170068","rs7177514","rs7183604","rs7359276","rs8042059","rs8042374","rs8042494","rs8043009","rs938682")
snp_to_gene <- getBM(attributes = c("refsnp_id", "allele", "chr_name",
                                    "chrom_start", "chrom_strand", "synonym_name",
                                    "ensembl_gene_stable_id", "consequence_type_tv"),
                     filters = "snp_filter",
                     values = rsids,
                     mart = ensembl)

gene_to_symbol <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                       filters ='ensembl_gene_id', 
                       values = snp_to_gene$ensembl_gene_stable_id, 
                       mart = genes)

mart_res <- merge(snp_to_gene, gene_to_symbol, by.x = "ensembl_gene_stable_id", 
                 by.y = "ensembl_gene_id", all.x = TRUE)

min_res <- mart_res %>% 
  dplyr::select("refsnp_id", "ensembl_gene_stable_id", 
         "hgnc_symbol", "consequence_type_tv") %>%
  distinct()

gwas <- merge(gwas, min_res, by.x = "SNP", by.y = "refsnp_id", 
              all.x = TRUE)
write_csv(gwas, "hits_overlap_annotated.csv")
