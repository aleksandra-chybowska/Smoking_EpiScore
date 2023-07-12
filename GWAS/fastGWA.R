#### Smoking as clean pack years

phenos = read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_17865_complete.csv")
head(phenos)
i = 8
tmp = data.frame(fid = phenos[,1],
	             iid = phenos[,1], 
	             phenp = phenos[,i])

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/transformed_pack_years")
write.table(tmp, file="phenotype.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(phenos[,c("Sample_Name", "Sample_Name")], file="samps.txt", sep='\t', col.names=F, row.names=F, quote=F)

./gcta64_dm  --mbfile mbfile_GS.txt \
--make-grm \
--sparse-cutoff 0.05 \
--threads 30 \
--out methid_grm \
--keep ./transformed_pack_years/samps.txt > GRM_gen.log


./gcta64_dm  --mbfile mbfile_GS.txt \
						  --grm-sparse methid_grm \
       					  --fastGWA-mlm \
       					  --maf 0.01 \
       					  --pheno ./transformed_pack_years/phenotype.phen \
       					  --threads 30 \
       					  --out ./transformed_pack_years/outfile


#### Smoking as epismoker
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/")

path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/"
samps = readRDS("../GS_methylation/GS20k/GS20k_Targets.rds")
epismk = read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/epismoker_20k.csv")

phenos = merge(samps, epismk, by.x="Sample_Sentrix_ID", by.y="SampleName")
table(is.na(phenos$smokingScore))
phenos = subset(phenos, !is.na(phenos$smokingScore))
i = 7
tmp = data.frame(fid = phenos[,2],
	             iid = phenos[,2], 
	             phenp = phenos[,i])

# phenos2 = read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_17865_complete.csv")
# head(phenos)
# i = 8
# tmp2 = data.frame(fid = phenos2[,1],
# 	             iid = phenos2[,1], 
# 	             phenp = phenos2[,i])
# identical(tmp$IID, tmp2$IID)
# IDs are identical, I can use the same grm for gwas

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/epismoker")
write.table(tmp, file="phenotype.phen", sep='\t', row.names=F, quote=F, col.names=F)
write.table(phenos[,c("Sample_Name", "Sample_Name")], file="samps.txt", sep='\t', col.names=F, row.names=F, quote=F)

./gcta64_dm  --mbfile mbfile_GS.txt \
						  --grm-sparse methid_grm \
       					  --fastGWA-mlm \
       					  --maf 0.01 \
       					  --pheno ./epismoker/phenotype.phen \
       					  --threads 30 \
       					  --out ./epismoker/outfile