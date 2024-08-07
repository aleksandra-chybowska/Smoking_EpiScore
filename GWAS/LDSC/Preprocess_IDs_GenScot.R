anno = read.csv('C:/Users/s1654019/Desktop/anno_full.csv', header=F)
epismoker = read.csv('C:/Users/s1654019/Desktop/outfile.fastGWA', sep='\t')

colnames(anno) = c("LOC", "rsID", "CHR", "POS", "A1", "A2", "FREQ", "ID")
saveRDS(anno, 'U:/Ola/Lab/Smoking/LDSC/anno.RDS', compress=F)


anno = anno[c("rsID", "ID")]
ds = merge(anno, epismoker, by.x = "ID", by.y = "SNP")
smoker = readRDS("C:/Users/s1654019/Desktop/EpiSmokEr_translated.outfile.fastGWA.RDS")
write.table(smoker, "C:/Users/s1654019/Desktop/EpiSmokEr_translated.outfile.fastGWA.tab", sep='\t', quote = F)


#======================

library(dplyr)
library(data.table)


#View(gwas_results)
setwd('/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS')
tpy = 'transformed_pack_years_overlap/outfile.fastGWA'

anno_full = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno_full, "V8", "V2")

tpy = read.table(tpy, header = T)
min_tpy = merge(tpy, anno, by.x="SNP", by.y="V8", all.x = T)
min_tpy$SNP = ifelse(is.na(min_tpy$V2), min_tpy$SNP, min_tpy$V2)
head(min_tpy)
fwrite(min_tpy, "GS_17105_pack_years_GWAS_translated_test.fastGWA", sep = "\t", row.names = F)