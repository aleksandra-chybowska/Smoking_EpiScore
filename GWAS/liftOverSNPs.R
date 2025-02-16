library(biomaRt)
library(data.table)
library(bigsnpr)
library(tidyverse)

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/EUR_stratified")
# read in Sauders

files = c("GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt",
  "GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt",
  "GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt",
  "GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt")

out <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/liftOver/"
liftover_path <- "/home/s1654019/micromamba/envs/openssl/bin/liftOver"

for (f in files) {
  file <- fread(f)
  colnames(file) = tolower(colnames(file))

  options(scipen = 999)

  bed = file %>% mutate(
    chromosome = chr,
    start = pos - 1,
    end = pos, 
    pasted = paste0(chromosome, "-", end)) %>%
    select(chromosome, start, end, pasted, rsid)

  fwrite(bed, 
    paste0(out, f, ".bed"), 
    sep="\t",  col.names = FALSE, quote = FALSE, row.names = FALSE)

  # Define input, output, and chain files
  input_bed <- paste0(out, f, ".bed")
  chain_file <- paste0(out, "hg38ToHg19.over.chain.gz")
  output_bed <- paste0(out, f, "hg19.bed")
  unmapped_bed <- paste0(out, f, "hg19.unmapped.bed")

  # Run liftOver using system2
  system2(liftover_path, 
          args = c(input_bed, chain_file, output_bed, unmapped_bed))

  output <- fread(paste0(out, f, "hg19.bed"))
  colnames(output) <- c("chromosome", "start", "stop", "pasted", "rsid")

  file$pasted = paste0(file$chr, "-", file$pos)
  merged <- merge(file, output, by="pasted")
  saveRDS(merged, paste0(out, f, ".merged.RDS"))

#  write_csv(merged, paste0(out, f, ".merged"))
}

for (f in files) {
  merged <- readRDS(paste0(out, f, ".merged.RDS"))
  # fix column naming and save
  target_cols <- c("chr", "pos", "rsid", "effect_allele", "other_allele",
    "af_1000g", "beta", "se", "p", "n")
  adapted <- c("chromosome", "stop", "rsid.x", "effect_allele", "other_allele",
    "af_1000g", "beta", "se", "p", "n")

  merged_hg19 = merged %>% 
    select(adapted) %>% 
    rename("rsid" = rsid.x, 
    "chr" = chromosome, 
    "pos" = stop
    )
  options(scipen = 999)
  fwrite(merged_hg19, 
    paste0(out, "complete/", f, ".hg19.txt"), 
    sep="\t", quote = FALSE, row.names = FALSE)
}


## BIOMART version - very slow 
rsid_to_coords <- function(mart, rsids) {
  # Query SNP coordinates in hg38 (GRCh38)
  snps <- getBM(
    attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end'),
    filters = 'snp_filter',
    values = rsids,
    mart = mart
  )

  return(snps)
}

# Connect to Ensembl and choose datasets
ensembl_hg19 <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh=37)
ensembl_hg38 <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh=38)


for (f in files) {
  file <- fread(f)
  rsids <- file$RSID

  snps_hg19 = list()
  for (i in seq(from=1, to=length(rsids), by=100001)) {
    batch_len = 100000
    if (length(rsids) < i + batch_len) {
      batch_len = length(rsids)
    }
    print(paste0("Start: ", i, ", stop:", i + batch_len))
    snps_hg19_batch <- rsid_to_coords(ensembl_hg19, rsids[i:batch_len])
    snps_hg19 <- c(snps_hg19, snps_hg19_batch)
  }

  print(f)
  print("Before translating from GRCh38, this was the length of summstats: ")
  print(dim(file))
  print("After translating:")
  print(dim(snps_hg19))
  print("After merge: ")

  translated = merge(file, snps_hg19, 
                      by.x = "RSID", by.y = "refsnp_id", all.x = TRUE)

  print("How many SNPs are missing?")
  print(table(is.na(translated$chrom_start)))
}
  # Query SNP coordinates in hg19 (GRCh37)

# Show the results
print(snps_hg19)
print(snps_hg38)