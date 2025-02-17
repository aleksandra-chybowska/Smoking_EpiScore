#library(dplyr)

#### Smoking as clean pack years
# phenos = read.csv("<cluster_home_dir>/Smoking/BayesRR/data/pack_years_17865_complete.csv")
# head(phenos)
# i = 8
# tmp = data.frame(fid = phenos[,1],
# 	             iid = phenos[,1], 
# 	             phenp = phenos[,i])

# setwd("<cluster_home_dir>/Smoking/GWAS/transformed_pack_years")
# write.table(tmp, file="phenotype.phen", sep='\t', row.names=F, quote=F, col.names=F)
# write.table(phenos[,c("Sample_Name", "Sample_Name")], file="samps.txt", sep='\t', col.names=F, row.names=F, quote=F)

# generate full (no keep flag) or project specific GRM (with keep flag)
# phenotype generation in a separate file

./gcta64_dm  --mbfile mbfile_GS.txt \
--make-grm \
--sparse-cutoff 0.05 \
--threads 30 \
--out methid_grm \
--keep ./transformed_pack_years/samps.txt > GRM_gen.log

./gcta64_dm  --mbfile mbfile_GS.txt \
					--grm-sparse methid_grm_full \
					--fastGWA-mlm \
					--maf 0.01 \
					--pheno ./transformed_pack_years_overlap/transformed_pack_years_overlap_17130.phen \
					--qcovar ./transformed_pack_years_overlap/qcovar_17130.phen \
					--covar ./transformed_pack_years_overlap/covar_17130.phen \
					--threads 30 \
					--out ./transformed_pack_years_overlap/outfile

./gcta64_dm  --mbfile mbfile_GS.txt \
					--grm-sparse methid_grm_full \
					--fastGWA-mlm \
					--maf 0.01 \
					--pheno ./epismoker_overlap/epismoker_overlap_17130.phen \
					--qcovar ./epismoker_overlap/qcovar_17130.phen \
					--covar ./epismoker_overlap/covar_17130.phen \
					--threads 30 \
					--out ./epismoker_overlap/outfile
					
					
./gcta64_dm  --mbfile mbfile_GS.txt \
					--grm-sparse methid_grm_full \
					--fastGWA-mlm \
					--maf 0.01 \
					--pheno ./raw/raw_17130.phen \
					--qcovar ./raw/qcovar_17130.phen \
					--covar ./raw/covar_17130.phen \
					--threads 30 \
					--out ./raw/outfile

./gcta64_dm  --mbfile mbfile_GS.txt \
					--grm-sparse methid_grm_full \
					--fastGWA-mlm \
					--maf 0.01 \
					--pheno ./inrt_overlap/inrt_overlap_17124.phen \
					--qcovar ./inrt_overlap/qcovar_17124.phen \
					--covar ./inrt_overlap/covar_17124.phen \
					--threads 30 \
					--out ./inrt_overlap/outfile

./gcta64_dm  --mbfile mbfile_GS.txt \
					--grm-sparse methid_grm_full \
					--fastGWA-mlm \
					--maf 0.01 \
					--pheno ./grimage_overlap/grimage_overlap_17130.phen \
					--qcovar ./grimage_overlap/qcovar_17130.phen \
					--covar ./grimage_overlap/covar_17130.phen \
					--threads 30 \
					--out ./grimage_overlap/outfile