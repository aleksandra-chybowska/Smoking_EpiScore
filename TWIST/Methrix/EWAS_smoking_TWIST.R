library("stringr")
library("foreach")
library("doParallel")

home = "/Cluster_Filespace/Marioni_Group/Ola/"
approach = "parallel"
filename = "pass15_40" 
#filename = "Covered_annotated"
#filename = "Covered_annotated_snp_rm"

#meth = readRDS(paste0(home, 'TWIST/Methrix/All_beds.RDS'))
#meth = readRDS(paste0(home, "TWIST/Methrix/processed.RDS")) # 40_10
meth = readRDS(paste0(home, 'TWIST/Methrix/Filtered/', filename, '.RDS')) # 989315 after filtering 831600
meth = na.omit(meth)
meth = t(meth)
# dim  > dim(meth)
# [1]      48 4720922 - ALL
# 		   48 3886109 - Covered
#		   48 2926003 - 40_10
# 

pheno = read.csv(paste0(home, "TWIST/Pheno/2022-07-28-Twist_Sample_Info.csv"))
pheno$smk = ifelse(pheno$ever_smoke == "Never", 0, 1)
table(pheno$smk)
pheno$methrix_id = paste0("X", pheno$Sample_Name)

#data checks
summary(pheno)
identical(rownames(meth), pheno$methrix_id) # FALSE
pheno = pheno[match(rownames(meth), pheno$methrix_id), ]
identical(rownames(meth), pheno$methrix_id) # TRUE 

table(is.na(pheno)) # OK 
table(as.factor(pheno$sex)) # OK
table(as.factor(pheno$smk)) # OK
table(is.na(meth)) # OK 

if (approach == "parallel") {

	logs <- paste0(home, "TWIST/Manual_EWAS/")
	cores <- detectCores()
	cl <- makeCluster(10, outfile = paste0(logs, "parallel_loops_", filename, "_snp_rm.txt"))
	registerDoParallel(cl)

	ewas <- foreach(
		i = 1:ncol(meth), 
		.combine = "rbind"
		) %dopar% {
			cpg = colnames(meth)[i]
			split = unlist(strsplit(cpg, "-", fixed = TRUE))
			chr = split[1]
			start = split[2]

			mod = summary(lm(meth[,i] ~ pheno$age + 
										as.factor(pheno$sex) + 
										as.factor(pheno$smk)))
			
			beta = mod$coefficients[4,1]
			se = mod$coefficients[4,2]
			tstat = mod$coefficients[4,3]
			p = mod$coefficients[4,4]

			return(data.frame(cpg, chr, start, beta, se, tstat, p))
	}

	saveRDS(ewas, paste0(home, "TWIST/Manual_EWAS/", filename, "_EWAS.RDS"))

} else {
	p = NA
	cpg = NA
	for(i in 1:ncol(meth)){
		p[i] = summary(lm(meth[,i] ~ pheno$age + as.factor(pheno$sex) + as.factor(pheno$smk)))$coefficients[4,4]
		cpg[i] = colnames(meth)[i]
	}

	idx = which(p < 1e-5)
	n = length(idx) # 97
	cpgs = cpg[idx]

	write.csv(idx, paste0(home, "TWIST/Manual_EWAS/idx.csv"), row.names=FALSE)
	write.csv(p[idx], paste0(home, "TWIST/Manual_EWAS/p.csv"), row.names=FALSE)
	write.csv(cpg[idx], paste0(home, "TWIST/Manual_EWAS/cpgs.csv"), row.names=FALSE)	
}