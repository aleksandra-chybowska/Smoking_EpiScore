library("qqman")

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/grimage_overlap")
res = read.table("outfile.fastGWA", header = T)

#this line converts each p-value to a chi-square value
# with 1 degree of freedom (typical for GWAS tests).
chi = qchisq(1-res$P,1)

# Lambda is calculated as the median of the chi-square values divided by 0.456,
# which is the expected median chi-square value under the null hypothesis 
# for a 1 degree of freedom test.
lambda = median(chi)/0.456
names(lambda) = "lambda"
write.csv(lambda, "lambda_tpy.csv", row.names = F)

min_ds = res[c("SNP", "CHR", "POS", "P")]

names(min_ds)<-c("SNP", "CHR", "BP", "P")
head(min_ds)

jpeg("qq.jpg")
qq(min_ds$P)
dev.off()

