library("qqman")

setwd("/Users/shirin/Documents/Edinburgh/Lab/Smoking/GWAS/transformed_pack_years_overlap")
res = read.table("outfile.fastGWA", header = T)
chi = qchisq(1-res$P,1)
lambda = median(chi)/0.456
names(lambda) = "lambda"
write.csv(lambda, "lambda_tpy.csv", row.names = F)

min_ds = res[c("SNP", "CHR", "POS", "P")]

names(min_ds)<-c("SNP", "CHR", "BP", "P")
head(min_ds)

jpeg("qq.jpg")
qq(min_ds$P)
dev.off()

