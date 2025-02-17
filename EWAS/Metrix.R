library(methrix)
library(BSgenome.Hsapiens.UCSC.hg38)

product = "unfiltered"
path = "<cluster_home_dir>/TWIST/Methrix/"

annotate = function(meth) {
  meth_mat_with_loci <- methrix::get_matrix(m = meth, type = "M", add_loci = TRUE, in_granges = TRUE)

  library(Repitools)
  anno = annoGR2DF(meth_mat_with_loci)
  rownames(anno)=paste0(anno$chr,"-",anno$start,"-",anno$end)
  anno = anno[-c(1:4)]
  head(colnames(anno))
  library(stringr)
  colnames(anno) = str_replace(colnames(anno), ".markdup.sorted_CpG.bedGraph", "")
  return(anno)
}

snp_rm = function(meth) {
  return(methrix::remove_snps(m = meth, keep = F))
}

bdg_files <- list.files(
  path = '<cluster_home_dir>/TWIST/Methylation',
  pattern = "*bedGraph",
  full.names = TRUE
)

bdg_files_test = bdg_files[1:4]

print(basename(bdg_files))

sample_anno <- data.frame(
  row.names = gsub(
    pattern = "\\.bedGraph\\$",
    replacement = "",
    x = basename(bdg_files)
  )
)

hg38_cpgs <- suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))

meth <- methrix::read_bedgraphs(
  files = bdg_files,
  pipeline = "MethylDackel",
  ref_cpgs = hg38_cpgs,
  coldata = sample_anno
)

# An object of class methrix
#    n_CpGs: 29,401,795
# n_samples: 48
#     is_h5: FALSE
# Reference: hg38

saveRDS(meth, paste0(path, 'All_beds.RDS'), compress = F)
saveRDS(annotate(meth), paste0(path, 'All_beds_annotated.RDS'), compress = F)

library(methrix)

# filter out sites of extremely low (untrustworthy) and high coverage (technical issue)
meth <- methrix::mask_methrix(meth, low_count = 2, high_quantile = 0.99) 

# filter to covered sites
meth <- methrix::remove_uncovered(meth)
# -Removed 11,153,323 [37.93%] uncovered loci of 29,401,795 sites
# -Finished in:  7.382s elapsed (6.006s cpu)

saveRDS(meth, paste0(path, 'Covered.RDS'), compress = F)
saveRDS(annotate(meth), paste0(path, 'Covered_annotated.RDS'), compress = F)
saveRDS(annotate(snp_rm(meth)), paste0(path, 'Covered_annotated_snp_rm.RDS'), compress = F)
# various coverage groups 
pass10 <- methrix::coverage_filter(m=meth, cov_thr = 10, min_samples = 48)
# -Retained 1,321,255 of 18,248,472 sites
# -Finished in:  3.958s elapsed (3.583s cpu)

pass10_44 <- methrix::coverage_filter(m=meth, cov_thr = 10, min_samples = 44)
# -Retained 2,806,718 of 18,248,472 sites
# -Finished in:  3.810s elapsed (3.614s cpu)

pass10_40 <- methrix::coverage_filter(m=meth, cov_thr = 10, min_samples = 40)
# -Retained 3,391,718 of 18,248,472 sites
# -Finished in:  3.823s elapsed (3.621s cpu)

pass10_36 <- methrix::coverage_filter(m=meth, cov_thr = 10, min_samples = 36)
# -Retained 3,769,095 of 18,248,472 sites
# -Finished in:  3.851s elapsed (3.637s cpu)

pass10_32 <- methrix::coverage_filter(m=meth, cov_thr = 10, min_samples = 32)
# -Retained 4,055,917 of 18,248,472 sites
# -Finished in:  3.850s elapsed (3.650s cpu)

pass15 <- methrix::coverage_filter(m=meth, cov_thr = 15, min_samples = 48)
# -Retained 578,617 of 18,248,472 sites
# -Finished in:  3.863s elapsed (3.649s cpu)
pass15_44 <- methrix::coverage_filter(m=meth, cov_thr = 15, min_samples = 44)
# -Retained 1,746,950 of 18,248,472 sites
# -Finished in:  3.860s elapsed (3.658s cpu)
pass15_40 <- methrix::coverage_filter(m=meth, cov_thr = 15, min_samples = 40)
# -Retained 2,332,188 of 18,248,472 sites
# -Finished in:  3.844s elapsed (3.631s cpu)
pass15_36 <- methrix::coverage_filter(m=meth, cov_thr = 15, min_samples = 36)
# -Retained 2,753,832 of 18,248,472 sites
# -Finished in:  3.838s elapsed (3.537s cpu)
pass15_32 <- methrix::coverage_filter(m=meth, cov_thr = 15, min_samples = 32)
# -Retained 3,096,514 of 18,248,472 sites
# -Finished in:  3.836s elapsed (3.595s cpu)

pass20 <- methrix::coverage_filter(m=meth, cov_thr = 20, min_samples = 48)
# Retained 216,596 of 18,248,472 sites
# -Finished in:  3.864s elapsed (3.658s cpu)
pass20_44 <- methrix::coverage_filter(m=meth, cov_thr = 20, min_samples = 44)
# -Retained 1,026,011 of 18,248,472 sites
# -Finished in:  3.845s elapsed (3.631s cpu)
pass20_40 <- methrix::coverage_filter(m=meth, cov_thr = 20, min_samples = 40)
# -Retained 1,517,055 of 18,248,472 sites
# -Finished in:  3.838s elapsed (3.613s cpu)
pass20_36 <- methrix::coverage_filter(m=meth, cov_thr = 20, min_samples = 36)
# -Retained 1,905,878 of 18,248,472 sites
# -Finished in:  3.835s elapsed (3.598s cpu)
pass20_32 <- methrix::coverage_filter(m=meth, cov_thr = 20, min_samples = 32)
# -Retained 2,245,234 of 18,248,472 sites
# -Finished in:  3.871s elapsed (3.657s cpu)

pass25 <- methrix::coverage_filter(m=meth, cov_thr = 25, min_samples = 48)
# -Retained 61,738 of 18,248,472 sites
# -Finished in:  3.876s elapsed (3.639s cpu)
pass25_44 <- methrix::coverage_filter(m=meth, cov_thr = 25, min_samples = 44)
# -Retained 570,157 of 18,248,472 sites
# -Finished in:  3.844s elapsed (3.619s cpu)
pass25_40 <- methrix::coverage_filter(m=meth, cov_thr = 25, min_samples = 40)
# -Retained 947,940 of 18,248,472 sites
# -Finished in:  3.829s elapsed (3.578s cpu)
pass25_36 <- methrix::coverage_filter(m=meth, cov_thr = 25, min_samples = 36)
# -Retained 1,266,416 of 18,248,472 sites
# -Finished in:  3.857s elapsed (3.572s cpu)
pass25_32 <- methrix::coverage_filter(m=meth, cov_thr = 25, min_samples = 32)
# -Retained 1,562,401 of 18,248,472 sites
# -Finished in:  3.849s elapsed (3.669s cpu)

pass30 <- methrix::coverage_filter(m=meth, cov_thr = 30, min_samples = 48)
# -Retained 11,099 of 18,248,472 sites
# -Finished in:  3.879s elapsed (3.663s cpu)
pass30_44 <- methrix::coverage_filter(m=meth, cov_thr = 30, min_samples = 44)
# -Retained 288,883 of 18,248,472 sites
# -Finished in:  4.098s elapsed (3.842s cpu)
pass30_40 <- methrix::coverage_filter(m=meth, cov_thr = 30, min_samples = 40)
# -Retained 570,131 of 18,248,472 sites
# -Finished in:  3.874s elapsed (3.623s cpu)
pass30_36 <- methrix::coverage_filter(m=meth, cov_thr = 30, min_samples = 36)
# -Retained 817,340 of 18,248,472 sites
# -Finished in:  3.847s elapsed (3.609s cpu)
pass30_32 <- methrix::coverage_filter(m=meth, cov_thr = 30, min_samples = 32)
# -Retained 1,056,618 of 18,248,472 sites
# -Finished in:  3.872s elapsed (3.598s cpu)

filtered = c(
  "pass10" = pass10, "pass10_44" = pass10_44, "pass10_40" = pass10_40, "pass10_36" = pass10_36, "pass10_32" = pass10_32,
  "pass15" = pass15, "pass15_44" = pass15_44, "pass15_40" = pass15_40, "pass15_36" = pass15_36, "pass15_32" = pass15_32,
  "pass20" = pass20, "pass20_44" = pass20_44, "pass20_40" = pass20_40, "pass20_36" = pass20_36, "pass20_32" = pass20_32,
  "pass25" = pass25, "pass25_44" = pass25_44, "pass25_40" = pass25_40, "pass25_36" = pass25_36, "pass25_32" = pass25_32,
  "pass30" = pass30, "pass30_44" = pass30_44, "pass30_40" = pass30_40, "pass30_36" = pass30_36, "pass30_32" = pass30_32
)

#sex chromosomes?

library(GenomicScores) 
library(MafDb.1Kgenomes.phase3.GRCh38)

# Cytosine (C) to thymine (T) polymorphisms are the most frequent transitions occurring in the human genome; often driven by the spontaneous deamination of a methylated cytosine at CpG dinucleotides to yield thymine. The appearance of mixed pyrimidines (C and T) at a single genomic location also parallels the outcome of chemistry used to detect differential DNA methylation. In the latter case, sodium bisulfite treatment converts unmethylated Cs to Ts, whereas methylated Cs remain unchanged. Thus, a common C to T polymorphism appearing specifically in one population could be misinterpreted as differential DNA methylation between individuals. Pinpointing sites across the array where this conflation may be occurring avoids erroneously calling differential methylation.

if (approach == "parallel") {

  library("foreach")
  library("doParallel")

	logs <- path
	cores <- detectCores()
	cl <- makeCluster(5, outfile = paste0(logs, "parallel_loops_eh.txt"))
	registerDoParallel(cl)

	foreach(
		file = filtered,
    name = names(filtered),
    .packages = c("methrix", "BSgenome.Hsapiens.UCSC.hg38")
		) %dopar% {
			file <- methrix::remove_snps(m = file, keep = F)
      file <- annotate(file)
      saveRDS(file, paste0(path, "Filtered/", name, ".RDS"), compress=F)
	}
}

if (F)
{


pass10_40 <- methrix::remove_snps(m = pass10_40, keep = F)

# 3,173,397
saveRDS(pass10_40, '<cluster_home_dir>/TWIST/Methrix/10_40_no_snps.RDS', compress=F)

plot = plot_coverage(pass10_40)
png("<cluster_home_dir>/TWIST/Methrix/Plots/coverage_10_40_no_snps.png", width = 1400, height=1400)
print(plot)
dev.off()

plot = plot_coverage(pass10_40, type="dens")
png("<cluster_home_dir>/TWIST/Methrix/Plots/coverage_10_40_no_snps_dens.png", width = 1400, height=1400)
print(plot)
dev.off()

library(data.frame)

meth_stats = get_stats(m = pass10_40)
fwrite(as.data.frame(meth_stats), file ="<cluster_home_dir>/TWIST/Methrix/pass10_40.csv")

meth_mat_with_loci <- methrix::get_matrix(m = pass10_40, type = "M", add_loci = TRUE, in_granges = TRUE)
saveRDS(meth_mat_with_loci, "<cluster_home_dir>/TWIST/Methrix/pass10_40_meth_granges.RDS", compress=FALSE)

coverage_mat_with_loci <- methrix::get_matrix(m = pass10_40, type = "C", add_loci = TRUE, in_granges = TRUE)
saveRDS(coverage_mat_with_loci, "<cluster_home_dir>/TWIST/Methrix/pass10_40_cov_granges.RDS", compress=FALSE)

library(Repitools)
anno = annoGR2DF(meth_mat_with_loci)
rownames(anno)=paste0(anno$chr,"-",anno$start,"-",anno$end)
anno = anno[-c(1:4)]
head(colnames(anno))
library(stringr)
colnames(anno) = str_replace(colnames(anno), ".markdup.sorted_CpG.bedGraph", "")

saveRDS(anno, "<cluster_home_dir>/TWIST/Methrix/processed.RDS", compress=FALSE)

#colnames(anno) = str_replace(colnames(anno), "X", "")
# > dim(anno)
# [1] 3173397      48
anno_no_nas = na.omit(anno)
# 2926003
} #SNP filtering 
# w1 after m2beta
#       201118840028_R03C01 201118840028_R04C01 201118840028_R05C01
# cg18478105          0.03111096          0.02958550          0.02693748
# cg09835024          0.05250425          0.05241946          0.03496545
# cg14361672          0.78746375          0.79019199          0.70569642
# cg01763666          0.72693236          0.69327342          0.72391360
# cg12950382          0.80162927          0.79261803          0.74779488
# cg02115394          0.10689030          0.11735072          0.12268549

# samples after QC

#                    X46261 X47853 X58855 X61901 X7394 X7467 X77847 X79037 X82551
# chr1-17406-17407     0.94   1.00   0.71   0.72  1.00  0.75   0.93   0.94   1.00
# chr1-17452-17453     0.97   0.78   0.68   0.96  0.97  0.80   1.00   1.00   0.46
# chr1-17478-17479     0.91   0.90   0.95   0.98  0.94  1.00   0.96   0.93   0.92
# chr1-17483-17484     0.93   0.97   0.95   0.98  0.94  0.95   0.93   0.95   1.00
# chr1-17562-17563     0.91   0.75   0.91   0.98  0.94  1.00   0.80   1.00   1.00
# chr1-135028-135029   1.00   1.00   0.95   0.94  0.91  0.92   0.93   1.00   1.00