library(stringr)

#read in methylation file
path = "/Users/shirin/Documents/Edinburgh/Lab/TWIST/Methylation/"
files = list.files(path, "*.bedGraph")
target = readRDS('/Users/shirin/Documents/Edinburgh/Lab/Methylation_test_sets/GS20k_Targets.rds')

full_w1 = nrow(read.csv(paste0(path, '../rows_w1.csv')))
full_w3 = nrow(read.csv(paste0(path, '../rows_w3.csv')))
full_w4 = nrow(read.csv(paste0(path, '../rows_w4.csv')))

stats = data.frame("sample" = character(), 
                   "wave"= character(), 
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
  
  file$pass_1  = ifelse(file$yes+file$no >= 1, T, F)
  file$pass_10 = ifelse(file$yes+file$no >= 10, T, F)
  file$pass_30 = ifelse(file$yes+file$no >= 30, T, F)
  file$pass_50 = ifelse(file$yes+file$no >= 50, T, F)
  
  wave = subset(target, Sample_Name == sample)$Set
  CpGs = ifelse(wave=="wave1", full_w1, ifelse(wave == "wave3", full_w3, full_w4))
  
  row = data.frame("sample" = sample, 
                   "wave" = wave,
                   "CpGs_850" = 860000,
                   "CpGs_wave" = CpGs,
                   "PASS1" = sum(file$pass_1),
                   "PASS10" = sum(file$pass_10),
                   "PASS30" = sum(file$pass_30),
                   "PASS50" = sum(file$pass_50)
                   )
  
  stats = rbind(stats, row)
}

write.csv(stats, paste0(path, "stats.csv"))

