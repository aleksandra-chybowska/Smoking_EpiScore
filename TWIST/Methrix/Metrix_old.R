library(methrix)
library(BSgenome.Hsapiens.UCSC.hg38)

bdg_files <- list.files(
  path = '/Cluster_Filespace/Marioni_Group/Ola/TWIST/Methylation',
  pattern = "*bedGraph",
  full.names = TRUE
)

bdg_files_test = bdg_files[1:4]

print(basename(bdg_files))

sample_anno <- data.frame(
  row.names = gsub(
    pattern = "\\.bedGraph\\$",
    replacement = "",
    x = basename(bdg_files_test)
  )
)

hg38_cpgs <- suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))

meth <- methrix::read_bedgraphs(
  files = bdg_files,
  pipeline = "MethylDackel",
  ref_cpgs = hg38_cpgs,
  coldata = sample_anno
)

meth <- methrix::read_bedgraphs(
  files = bdg_files[1:4],
  pipeline = "MethylDackel",
  ref_cpgs = hg38_cpgs,
  coldata = sample_anno,
  stranded = FALSE,
  zero_based = TRUE,
)
meth <- methrix::read_bedgraphs(
  files = bdg_files[1:4],
  ref_cpgs = hg38_cpgs,
  chr_idx = 1,
  start_idx = 2,
  end_idx = 3,
  M_idx = 5,
  U_idx = 6, 
  coldata = sample_anno,
  stranded = FALSE,
  zero_based = TRUE
)
