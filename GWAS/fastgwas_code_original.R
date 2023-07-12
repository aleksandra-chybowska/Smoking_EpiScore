# Prepare phenotypes in R
phenos = # phenotype file with id, phenotype columns
head(phenos)

tmp = data.frame(fid = phenos[,1],
	             iid = phenos[,1], 
	             phenp = phenos[,i])
write.table(tmp, file="phenotype.phen"), sep='\t', row.names=F, quote=F, col.names=F)


write.table(phenos[,c("ID", "ID")], file="samps.txt", sep='\t', col.names=F, row.names=F, quote=F)


# # Make sparse GRM
 ./gcta64_dm  --mbfile gs_mbfile.txt --make-grm --sparse-cutoff 0.05 --threads 10 --out methid_grm --keep samps.txt


# Run fastgwa, fitting sparse GRM
# mkdir fastGWA
./gcta64_dm  --mbfile gs_mbfile.txt \
						  --grm-sparse methid_grm \
       					  --fastGWA-mlm \
       					  --maf 0.01 \
       					  --pheno ./phenotype.phen \
       					  --threads 30 \
       					  --out ./outfile

