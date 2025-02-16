# Sex-agnostic
datadir <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/genetics_W3/"
filename <- "pack_years_4375_wave3_genetic_age_sex.csv"


pack_years <- read.csv(paste0(datadir, '../' ,filename))
cpgs <- fread(paste0(datadir, "../methylation_W3_genetic/GS20k_allchrom.csv"))
cpgs <- as.data.frame(cpgs)

# add Age and Sex from pheno
agesex <- data.frame(age = log(pack_years$age), sex = pack_years$sex)
rownames(agesex) <- agesex$Sample_Name
agesex$sex <- ifelse(agesex$sex == "F", 1, 0)
agesex <- scale(agesex)
names <- as.data.frame(colnames(agesex))
colnames(names) <- c("Marker")

agesex <- t(agesex)
data <- rbind(cpgs, agesex)

fwrite(x = data, paste0(datadir, "../multiomics/pack_years_W3_4375_agesex_meth_only_processed.csv"), 
	sep = ",", row.names = F, col.names = F, quote = F) 
rm(cpgs, data)
gc()

# Fuse SNPs and CpGs
cpgs_list <- read.csv(paste0(datadir, "../methylation_W3_genetic/GS20k_allchrom_cpg_list.txt"))
markers <- rbind(cpgs_list, names)

fwrite(x = as.matrix(markers), paste0(datadir, "../multiomics/pack_years_snp_list_meth_agesex.txt"), 
	sep = ",", row.names = F, col.names = T, quote = F) 
