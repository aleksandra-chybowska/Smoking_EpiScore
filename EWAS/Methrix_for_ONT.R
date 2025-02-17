#Load library
library(methrix)

#Genome of your preference to work with
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38")) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
}
if (!requireNamespace("Repitools")) {
  BiocManager::install("Repitools")
}
if (!requireNamespace("MafDb.1Kgenomes.phase3.GRCh38")) {
  BiocManager::install("MafDb.1Kgenomes.phase3.GRCh38")
}

if (!requireNamespace("GenomicScores")) {
  BiocManager::install("GenomicScores")
}

library("Repitools")
library("BSgenome.Hsapiens.UCSC.hg38")

setwd("<cluster_home_dir>/ONT/final_bedgraphs")

bdg_files <- list.files(
  path = ".",
  pattern = "*bedmethyl.gz",
  full.names = TRUE
)
# test
# bdg_files = c("./2947.wf_mods.no_WS.bedmethyl")

hg38_cpgs <- suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))

meth <- methrix::read_bedgraphs(
  files = bdg_files,
  ref_cpgs = hg38_cpgs,
  chr_idx = 1,
  start_idx = 2,
  strand_idx = 6,
  end_idx = 3,
  cov_idx = 10,
  beta_idx = 11,
  zero_based = TRUE,
  collapse_strands = FALSE,
  stranded = FALSE
)

saveRDS(meth, "../methrix/48_bedgraphs_unprocessed.RDS", compress = F)
meth = readRDS("../methrix/48_bedgraphs_unprocessed.RDS")

# F2RL3 = c("chr19-16889741-16889742",
#         "chr19-16889774-16889775",
#         "chr19-16889785-16889786")
test_F2RL3 <- function(meth) {
  cov_mat <- methrix::get_matrix(m = meth, type = "C", 
    add_loci = TRUE, in_granges = TRUE)
  meth_mat <- methrix::get_matrix(m = meth, type = "M", 
    add_loci = TRUE, in_granges = TRUE)

  queries <- lapply(list(c(16889741, 16889742), 
                        c(16889774, 16889775), 
                        c(16889785, 16889786)), 
                    function(coords) {
                      GRanges(seqnames="chr19",
                      ranges=IRanges(start=coords[1], end=coords[2]),
                      strand="*")})

  # Apply subsetByOverlaps to each query
  result_meth <- lapply(queries, function(q) 
    subsetByOverlaps(meth_mat, q))
  result_cov <- lapply(queries, function(q) 
    subsetByOverlaps(cov_mat, q))
  return(c(result_meth, result_cov))
}


#An object of class methrix
#    n_CpGs: 29,401,795
# n_samples: 30
#     is_h5: FALSE
# Reference: hg38

# filter out sites of extremely low (untrustworthy) and high coverage (technical issue)
meth2 <- methrix::mask_methrix(meth, low_count = 2, high_quantile = 0.99)

# filter to covered sites
meth2 <- methrix::remove_uncovered(meth2)
# Removed 412,393 [1.4%] uncovered loci of 29,401,795 sites
# -Finished in:  11.0s elapsed (8.533s cpu)


saveRDS(meth2, "../methrix/48_bedgraphs_covered.RDS", compress = F)
meth2 = readRDS("../methrix/48_bedgraphs_covered.RDS")

annotate = function(meth) {
  meth_mat_with_loci <- methrix::get_matrix(m = meth, type = "M", add_loci = TRUE, in_granges = TRUE)

  library(Repitools)
  anno = annoGR2DF(meth_mat_with_loci)
  rownames(anno)=paste0(anno$chr,"-",anno$start,"-",anno$end)
  anno = anno[-c(1:4)]
  head(colnames(anno))
  library(stringr)
  colnames(anno) = str_replace(colnames(anno), "collapsed_", "")
  return(anno)
}

snp_rm = function(meth) {
  return(methrix::remove_snps(m = meth, keep = F))
}

pass1_48 <- methrix::coverage_filter(m = meth2, cov_thr = 1, min_samples = 48)
# -Retained 13,734,218 of 28,989,402 sites
# -Finished in:  9.666s elapsed (8.204s cpu)
pass1_44 <- methrix::coverage_filter(m = meth2, cov_thr = 1, min_samples = 44)
# -Retained 26,661,741 of 28,989,402 sites
# -Finished in:  11.1s elapsed (8.763s cpu)
pass1_40 <- methrix::coverage_filter(m = meth2, cov_thr = 1, min_samples = 40)
# -Retained 27,352,245 of 28,989,402 sites
# -Finished in:  12.4s elapsed (9.679s cpu)
pass1_36 <- methrix::coverage_filter(m = meth2, cov_thr = 1, min_samples = 36)
# -Retained 27,660,684 of 28,989,402 sites
# -Finished in:  14.1s elapsed (10.1s cpu)
pass1_32 <- methrix::coverage_filter(m = meth2, cov_thr = 1, min_samples = 32)
# -Retained 27,874,117 of 28,989,402 sites
# -Finished in:  12.8s elapsed (9.409s cpu)



pass5_48 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 48)
# -Retained 148,623 of 28,989,402 sites
# -Finished in:  9.815s elapsed (8.440s cpu)

pass5_44 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 44)
# -Retained 9,760,346 of 28,989,402 sites
# -Finished in:  12.7s elapsed (9.577s cpu)

pass5_40 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 40)
# -Retained 21,167,712 of 28,989,402 sites
# -Finished in:  9.632s elapsed (8.116s cpu)

pass5_36 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 36)
# -Retained 24,879,584 of 28,989,402 sites
# -Finished in:  14.3s elapsed (10.4s cpu)

pass5_32 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 32)
# -Retained 26,399,654 of 28,989,402 sites
# -Finished in:  12.4s elapsed (9.341s cpu)

pass5_28 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 28)
# -Retained 27,186,136 of 28,989,402 sites
# -Finished in:  10.2s elapsed (8.787s cpu)

pass5_25 <- methrix::coverage_filter(m = meth2, cov_thr = 5, min_samples = 25)
# -Retained 27,510,631 of 28,989,402 sites
# -Finished in:  9.305s elapsed (8.077s cpu)
# ONT recommend this, maybe with even greater min_samples


# closest to TWIST - approx 1, 000 000 sites
pass7_48 <- methrix::coverage_filter(m = meth2, cov_thr = 7, min_samples = 48)
# -Retained 57 of 28,989,402 sites
# -Finished in:  13.0s elapsed (9.772s cpu)

pass7_44 <- methrix::coverage_filter(m = meth2, cov_thr = 7, min_samples = 44)
# -Retained 163,938 of 28,989,402 sites
# -Finished in:  12.7s elapsed (9.633s cpu)

pass7_40 <- methrix::coverage_filter(m = meth2, cov_thr = 7, min_samples = 40)
# -Retained 3,642,242 of 28,989,402 sites
# -Finished in:  9.304s elapsed (8.022s cpu)
# This is really closest to TWIST

pass7_36 <- methrix::coverage_filter(m = meth2, cov_thr = 7, min_samples = 36)
# -Retained 12,790,554 of 28,989,402 sites
# -Finished in:  12.5s elapsed (9.421s cpu)

pass7_32 <- methrix::coverage_filter(m = meth2, cov_thr = 7, min_samples = 32)
# -Retained 20,182,538 of 28,989,402 sites
# -Finished in:  12.9s elapsed (9.489s cpu)


pass10_48 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 48)
# -Retained 0 of 28,989,402 sites
# -Finished in:  12.5s elapsed (9.648s cpu)

pass10_44 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 44)
# -Retained 0 of 28,989,402 sites
# -Finished in:  14.4s elapsed (10.4s cpu)


pass10_40 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 40)
# -Retained 471 of 28,989,402 sites
# -Finished in:  9.667s elapsed (8.251s cpu)


pass10_36 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 36)
# -Retained 19,417 of 28,989,402 sites
# -Finished in:  12.4s elapsed (9.568s cpu)

pass10_32 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 32)


pass10_25 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 25)
# -Retained 8,505,191 of 28,989,402 sites
# -Finished in:  9.522s elapsed (8.146s cpu)

pass10_10 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 10)
# -Retained 26,845,960 of 28,989,402 sites
# -Finished in:  9.532s elapsed (8.126s cpu)

pass10_20 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 20)
# -Retained 18,142,315 of 28,989,402 sites
# -Finished in:  9.198s elapsed (7.934s cpu)

pass10_30 <- methrix::coverage_filter(m = meth2, cov_thr = 10, min_samples = 30)
# -Retained 1,358,510 of 28,989,402 sites
# -Finished in:  9.269s elapsed (7.994s cpu)

# Three interesting files:
pass5_40 <- methrix::remove_snps(m = pass5_40, keep = F)
pass1_40 <- methrix::remove_snps(m = pass1_40, keep = F)
pass5_36 <- methrix::remove_snps(m = pass5_36, keep = F)

saveRDS(pass5_40, "../methrix/5_40_raw_no_snps.RDS", compress = F)
saveRDS(pass1_40, "../methrix/1_40_raw_no_snps.RDS", compress = F)
saveRDS(pass5_36, "../methrix/5_36_raw_no_snps.RDS", compress = F)

file_annotated_pass5_40 = annotate(pass5_40)
file_annotated_pass1_40 = annotate(pass1_40)
file_annotated_pass5_36 = annotate(pass5_36)

# check whats inside file before running it!
saveRDS(file_annotated_pass5_40, "../methrix/5_40_after_methrix.RDS", compress = F)
dim(file_annotated_pass5_40)
saveRDS(file_annotated_pass1_40, "../methrix/1_40_after_methrix.RDS", compress = F)
dim(file_annotated_pass1_40)
saveRDS(file_annotated_pass5_36, "../methrix/5_36_after_methrix.RDS", compress = F)
dim(file_annotated_pass5_36)
