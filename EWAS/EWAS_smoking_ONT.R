library("stringr")
library("foreach")
library("doParallel")

setwd("<cluster_home_dir>/ONT/methrix")
home = "<cluster_home_dir>/"
approach = "normal"
filename = "5_40" 

meth = readRDS(paste0(filename, "_after_methrix.RDS"))
# > dim(meth)
# [1] 3642242      48 # 7_40

meth = na.omit(meth)
# [1] 2468919      48
meth = t(meth)
# [1]      48 2468919

pheno = read.csv(paste0(home, "TWIST/Pheno/2024-11-14-Twist_Sample_Info_Cleaned.csv"))
pheno$smk = ifelse(pheno$ever_smoke == "Never", 0, 1)
table(pheno$smk)
# pheno$methrix_id = paste0("X", pheno$Sample_Name)
pheno$methrix_id = paste0(pheno$Sample_Name)


#data checks
summary(pheno)
identical(rownames(meth), pheno$methrix_id) # FALSE
#remove additional rows
meth = subset(meth, rownames(meth) %in% pheno$methrix_id)
pheno = pheno[match(rownames(meth), pheno$methrix_id), ]
identical(rownames(meth), pheno$methrix_id) # TRUE 

table(is.na(pheno)) # OK 
table(as.factor(pheno$sex)) # OK
table(as.factor(pheno$smk)) # OK
table(is.na(meth)) # OK 

saveRDS(meth, paste0(filename, "_before_EWAS.RDS"), compress = F)
saveRDS(pheno, paste0(filename, "_before_EWAS_pheno.RDS"), compress = F)

### EWAS

meth = readRDS(paste0(filename, "_before_EWAS.RDS"))
pheno = readRDS(paste0(filename, "_before_EWAS_pheno.RDS"))

if (approach == "parallel") {

    logs <- paste0(home, "TWIST/Manual_EWAS/ont/")
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

    saveRDS(ewas, paste0(home, "TWIST/Manual_EWAS/ont/", filename, "_EWAS.RDS"))

} else {
    p = NA
    cpg = NA
    beta = NA
    se = NA
    tstat = NA

    for(i in 1:ncol(meth)){
        mod = summary(lm(meth[,i] ~ pheno$age +
                                    as.factor(pheno$sex) +
                                    as.factor(pheno$smk)))
        beta[i] = mod$coefficients[4,1]
        se[i] = mod$coefficients[4,2]
        tstat[i] = mod$coefficients[4,3]
        p[i] = mod$coefficients[4,4]
        cpg[i] = colnames(meth)[i]
    }

    split_data1 <- strsplit(cpg, "-")
    df1 <- do.call(rbind, split_data1)

    results = data.frame(cpg = cpg, chr = df1[,1], 
                        start=df1[,2], beta = beta, 
                        se = se, tstat = tstat, p = p)

    # idx = which(p < 0.05)
    # n = length(idx)
    # cpgs = cpg[idx]
    # setwd("<cluster_home_dir>/ONT/EWAS")

    # split_data <- strsplit(cpgs, "-")
    # df <- do.call(rbind, split_data)


    # # only marginally significant
    # results = data.frame(cpg = cpgs, chr = df[,1], 
    #                     start=df[,2], beta = beta[idx], 
    #                     se = se[idx], tstat = tstat[idx], p = p[idx])

    # all


    write.csv(results, paste0(filename, "46_full_EWAS_ONT.csv"), row.names=FALSE)
}
