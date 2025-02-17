library(stringr)

home = "<cluster_home_dir>/"

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

## Wave 3
meth_w3 = readRDS("<filespace_marioni_group_dir>/GS/GS_methylation/wave3_mvals.rds")
meth_w3[which(is.nan(meth_w3))] <- NA
meth_w3[which(is.infinite(meth_w3))] <- NA

# Mean impute missings
w3_nas = apply(meth_w3,1,function(x){which(is.na(x))})
nas = which(unlist(lapply(w3_nas, length))>0)
for(i in names(nas)){
  imp = mean(meth_w3[i,], na.rm=T)
  ind = which(is.na(meth_w3[i,]))
  meth_w3[i,ind] = rep(imp, length(ind))
}

which(is.na(meth_w3))
# numeric(0)
which(is.infinite(meth_w3))
# numeric(0)

# are the sites from TWIST significant?
# cpgs_twist = read.csv(paste0(home, "TWIST/Manual_EWAS/cpgs.csv"))
# this makes no sense since cpgs_twist[i] is a GRange

# list_w3 = read.table("<filespace_marioni_group_dir>/Daniel/GS_Wave3_EWAS/w3_probes_to_keep.txt")
# meth_w3 = meth_w3[which(row.names(meth_w3) %in% list_w3$V1),] 
meth = t(meth_w3)
meth = m2beta(meth)
rm(meth_w3)
gc()

pheno = read.csv(paste0(home, "TWIST/Pheno/2024-11-14-Twist_Sample_Info_Cleaned.csv"))
pheno$smk = ifelse(pheno$ever_smoke == "Never", 0, 1)
table(pheno$smk)
pheno$methrix_id = paste0("X", pheno$Sample_Name)

#data checks
summary(pheno)
# Now subset the DNAm by phenotype file IDs 
meth = meth[which(rownames(meth) %in% pheno$Sample_Sentrix_ID), ]
identical(rownames(meth), pheno$Sample_Sentrix_ID) # FALSE
dim(pheno)
pheno = pheno[match(rownames(meth), pheno$Sample_Sentrix_ID), ]
identical(rownames(meth), pheno$Sample_Sentrix_ID) # TRUE 

table(is.na(pheno)) # OK 
table(as.factor(pheno$sex)) # OK
table(as.factor(pheno$smk)) # OK
table(is.na(meth)) # OK 

# TODO: https://stackoverflow.com/questions/42400389/how-to-output-two-vectors-that-are-iteratively-filled-using-foreach
p = NA
cpg = NA

data <- data.frame(cpg = character(),
                  chromosome = numeric(),
                  position = numeric(),
                  beta = numeric(),
                  se = numeric(),
                  tstat = numeric(),
                  pval = numeric(),
                  study = character()
                )

for(i in 1:ncol(meth)){
  summary = summary(lm(meth[,i] ~ pheno$age + as.factor(pheno$sex) + as.factor(pheno$smk)))

	p[i] = summary(lm(meth[,i] ~ pheno$age + as.factor(pheno$sex) + as.factor(pheno$smk)))$coefficients[4,4]
	cpg[i] = colnames(meth)[i]

  test <- data.frame("cpg" = colnames(meth)[i],
                  "chromosome" = 0,
                  "position" = 0,
                  "beta" = summary$coefficients[4,1],
                  "se" = summary$coefficients[4,2],
                  "tstat" = summary$coefficients[4,3],
                  "pval" = summary$coefficients[4,4],
                  "study" = "W3")

  data <- rbind(test, data)
}

write.csv(data, paste0(home, "TWIST/Manual_EWAS/W3_46.csv"), row.names=FALSE)


idx = which(p < 1e-5)
n = length(idx) # 97
cpgs = cpg[idx]

write.csv(idx, paste0(home, "TWIST/Manual_EWAS/W3_idx.csv"), row.names=FALSE)
write.csv(p[idx], paste0(home, "TWIST/Manual_EWAS/W3_p.csv"), row.names=FALSE)
write.csv(cpg[idx], paste0(home, "TWIST/Manual_EWAS/W3_cpgs.csv"), row.names=FALSE)
