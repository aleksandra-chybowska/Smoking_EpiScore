library(dplyr)
cat = read.csv('<local_home_dir>/Downloads/smoking(4).tsv', sep='\t')
pip80 = read.csv('<local_home_dir>/Desktop/PIP_80.csv')
cat = cat %>% arrange(desc(n))
cat = subset(cat, !duplicated(cat$cpg))
cat = cat[c("cpg", "n", "author", "pmid", "p")]
cat = merge(pip80, cat, by.x="CpG.site", by.y="cpg", all.x=T)
cat = cat %>% arrange(desc(PIP))
write.csv(cat, '<local_home_dir>/Desktop/cat.csv', row.names = F)

sig = subset(cat, PIP>0.95 & p<1e-7)
