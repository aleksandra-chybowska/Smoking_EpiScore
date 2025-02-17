library(stringr)

#read in methylation files (bedGraph)
path = "<cluster_home_dir>/TWIST/Methylation/"
convert = "<cluster_home_dir>/TWIST/RDS_methylation/"
files = list.files(path, "*.bedGraph")

stats = data.frame("sample" = character(),
                   "wave" = character(),
                   "CpGs_850" = numeric(),
                   "CpGs_wave" = numeric(),
                   "PASS1" = numeric(),
                   "PASS10" = numeric(),
                   "PASS30" = numeric(),
                   "PASS50" = numeric()
                   )

for (hook in files) {
  print(hook)
  sample = unlist(str_split(hook, "\\."))[1]
  file = read.csv(file=paste0(path, hook), header=F, sep="\t")
  file = file[-1,]
  colnames(file) = c("chr", "start", "stop", "meth", "yes", "no")
  saveRDS(file, paste0(convert, sample, ".rds"), compress=F)

  file$pass_1  = ifelse(file$yes+file$no >= 1, T, F)
  file$pass_10 = ifelse(file$yes+file$no >= 10, T, F)
  file$pass_30 = ifelse(file$yes+file$no >= 30, T, F)
  file$pass_50 = ifelse(file$yes+file$no >= 50, T, F)

  row = data.frame("sample" = sample,
                   "wave" = "wave3",
                   "CpGs_850" = 860000,
                   "CpGs_wave" = "773860",
                   "PASS1" = sum(file$pass_1),
                   "PASS10" = sum(file$pass_10),
                   "PASS30" = sum(file$pass_30),
                   "PASS50" = sum(file$pass_50)
                   )

  stats = rbind(stats, row)
}

write.csv(stats, paste0(path, "stats.csv"))
