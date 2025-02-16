#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena (based on Rob's and Daniel's scripts)

library("data.table")
library("foreach")
library("doParallel")

path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/results/runs/white_blood_cells_complete/"
methylation = "methylation_complete"
priors = 4
# (1) Check each chain for convergence (rowSums of sigmas should be ~1)
###################################################################

setwd(paste0(path, "sigma/"))


pdf("check_convergence_sigma.pdf")
file = list.files("./", pattern = "*.csv")[1]
tmp = read.csv(file)
plot(rowSums(tmp), ylab="Sigma", main="Convergence")
dev.off()

# (2) Calculate AC variance explained by epigenetic probes
###################################################################

setwd(path)
loop = list.files("comp/", pattern = ".csv")
names = read.csv(paste0("../../../data/", methylation, "/GS20k_allchrom_cpg_list.txt"), header=T)
#names = read.csv("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/multiomics/pack_years_snp_list_meth_agesex.txt")
for(file in loop){
  output <- matrix(nrow = 1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Model"
  output$Model <- gsub(".csv*", "", file)
  
  ### (1) Mean variance explained by epigenetic probes and credible intervals
  sigma <- read.csv(paste0("sigma/", file))  
  output$Epigenetic_Mean_Variance_Explained <- mean(sigma[,2]/rowSums(sigma[,1:2]), na.rm = T)
  output$Epigenetic_Low_CI <- quantile(sigma[,2]/rowSums(sigma[,1:2]), na.rm = T, prob = 0.025)
  output$Epigenetic_High_CI <- quantile(sigma[,2]/rowSums(sigma[,1:2]), na.rm = T, prob = 0.975)
  
  ### (2) Calculate the proportion of variance that is attributable to small, medium and large effects - 1,2,3
  betas <- fread(paste0("beta/", file))
  betas <- as.data.frame(betas)
  comp <- fread(paste0("comp/", file))
  comp <- as.data.frame(comp)
  names(comp) <- names$Marker
  names(betas) <- names$Marker
  list <- apply(comp, 1, function(x) which(!x %in% c(1,2,3,4)))

  len = nrow(comp)

  x <- as.matrix(0, ncol = 1, nrow = len)
  x <- as.data.frame(x) # this transposes data
  for (i in 1:len) { 
    x[[i]] <- length(list[[i]]) == ncol(comp)
  }

  if (length(which(x[1,] %in% "TRUE")) > 0) { 
    comp <- comp[-which(x %in% "TRUE"),]
  } else { 
    comp <- comp
  }

  # Probes of small effect
  t <- vector() 
  list <- apply(comp, 1, function(x) which(x %in% 1))

# E.G.
# [[4]]
#   [1]    832   3012   3315   5513  13029  13657  13667  16814  18826  20068
#  [11]  20203  29743  30278  31712  40412  41331  44468  45668  48477  51223

  for (i in 1:len) { 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true <- which(t %in% "TRUE")
  ind_false <- which(t %in% "FALSE")
  list_true <- list[ind_true]
  list_false <- list[ind_false] 
  n <- length(list_true) 
  m1_1_true <- matrix(0, ncol = 1, nrow = n)
  m1_1_true <- as.data.frame(m1_1_true) 
  m1_1_true$ind <- ind_true
  x <- vector()
  for (j in m1_1_true$ind) { 
    x[j] <- sum((betas[j, list[[j]]])^2) 
  } 
  m1 <- as.data.frame(x) 
  m1$x[is.na(m1$x)] <- 0 
  names(m1) <- "Variance_Small_Effects" 
  
  # Probes of medium effect
  t <- vector() 
  list <- apply(comp, 1, function(x) which(x %in% 2)) 
  for (i in 1:len) { 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true <- which(t %in% "TRUE")
  ind_false <- which(t %in% "FALSE")
  list_true <- list[ind_true]
  list_false <- list[ind_false] 
  n <- length(list_true) # 997
  m2_true <- matrix(0, ncol = 1, nrow = n)
  m2_true <- as.data.frame(m2_true) 
  m2_true$ind <- ind_true
  x <- vector()
  for (j in m2_true$ind) { 
    x[j] <- sum((betas[j, list[[j]]])^2) 
  } 
  m2 <- as.data.frame(x) 
  m2$x[is.na(m2$x)] <- 0 
  names(m2) <- "Variance_Medium_Effects"
  
  # Probes of large effect
  t <- vector() 
  list <- apply(comp, 1, function(x) which(x %in% 3)) 
  for (i in 1:len){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true <- which(t %in% "TRUE")
  ind_false <- which(t %in% "FALSE")
  list_true <- list[ind_true]
  list_false <- list[ind_false] 
  n <- length(list_true) 
  m3_true <- matrix(0, ncol = 1, nrow = n)
  m3_true <- as.data.frame(m3_true) 
  m3_true$ind <- ind_true
  x <- vector()
  for (j in m3_true$ind) { 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m3 <- as.data.frame(x) 
  m3$x[is.na(m3$x)] <- 0 
  names(m3) <- "Variance_Large_Effects"

  if (priors >= 4) {
    # Probes of extra large effect
    t <- vector() 
    list <- apply(comp, 1, function(x) which(x %in% 4)) 

    for (i in 1:len){ 
      t[[i]] <- length(list[[i]]) > 0 
    } 
    ind_true <- which(t %in% "TRUE")
    ind_false <- which(t %in% "FALSE")
    list_true <- list[ind_true]
    list_false <- list[ind_false] 
    n <- length(list_true) 
    m4_true <- matrix(0, ncol = 1, nrow = n)
    m4_true <- as.data.frame(m4_true) 
    m4_true$ind <- ind_true
    x <- vector()
    for (j in m4_true$ind) { 
      x[j] <- sum((betas[j,list[[j]]])^2) 
    } 
    m4 <- as.data.frame(x) 
    m4$x[is.na(m4$x)] <- 0 
    names(m4) <- "Variance_Extra_Large_Effects"
  
    m4$num <- row.names(m4)
  }
  
  if (priors >= 5) {
    # Probes of extra large effect
    t <- vector() 
    list <- apply(comp, 1, function(x) which(x %in% 5)) 

    for (i in 1:len){ 
      t[[i]] <- length(list[[i]]) > 0 
    } 
    ind_true <- which(t %in% "TRUE")
    ind_false <- which(t %in% "FALSE")
    list_true <- list[ind_true]
    list_false <- list[ind_false] 
    n <- length(list_true) 
    m5_true <- matrix(0, ncol = 1, nrow = n)
    m5_true <- as.data.frame(m5_true) 
    m5_true$ind <- ind_true
    x <- vector()
    for (j in m5_true$ind) { 
      x[j] <- sum((betas[j,list[[j]]])^2) 
    } 
    m5 <- as.data.frame(x) 
    m5$x[is.na(m5$x)] <- 0 
    names(m5) <- "Variance_Extra_Large_Effects"
  
    m5$num <- row.names(m5)
  }

  if (priors == 3) {
   # Fuse
    m1$num <- row.names(m1) 
    m2$num <- row.names(m2) 
    m3$num <- row.names(m3) 
    all <- merge(m1, m2, by = "num", all = T) 
    var <- merge(all, m3, by = "num", all = T) 
    var[is.na(var)] <- 0 
    var$num <- NULL
    var$Total_Variance <- var[,1] + var[,2] + var[,3]
    var$Proportion_Small_Effects <- var[,1]/var[,4]
    var$Proportion_Medium_Effects <- var[,2]/var[,4]
    var$Proportion_Large_Effects <- var[,3]/var[,4]
    output$Proportion_Small_Effects <- mean(var$Proportion_Small_Effects) 
    output$Proportion_Medium_Effects <- mean(var$Proportion_Medium_Effects) 
    output$Proportion_Large_Effects <- mean(var$Proportion_Large_Effects) 
  }
  if(priors == 4) {
    # Fuse
    m1$num <- row.names(m1) 
    m2$num <- row.names(m2) 
    m3$num <- row.names(m3) 
    m4$num <- row.names(m4) 
    all <- merge(m1, m2, by = "num", all = T) 
    var <- merge(all, m3, by = "num", all = T) 
    var <- merge(var, m4, by = "num", all = T) 
    var[is.na(var)] <- 0 
    var$num <- NULL
    var$Total_Variance <- var[,1] + var[,2] + var[,3] + var[,4]
    var$Proportion_Small_Effects <- var[,1]/var[,5]
    var$Proportion_Medium_Effects <- var[,2]/var[,5]
    var$Proportion_Large_Effects <- var[,3]/var[,5]
    var$Proportion_Extra_Large_Effects <- var[,4]/var[,5]
    output$Proportion_Small_Effects <- mean(var$Proportion_Small_Effects) 
    output$Proportion_Medium_Effects <- mean(var$Proportion_Medium_Effects) 
    output$Proportion_Large_Effects <- mean(var$Proportion_Large_Effects) 
    output$Proportion_Extra_Large_Effects <- mean(var$Proportion_Extra_Large_Effects)
  }
  if(priors == 5) {
    # Fuse
    m1$num <- row.names(m1) 
    m2$num <- row.names(m2) 
    m3$num <- row.names(m3) 
    m4$num <- row.names(m4)
    m5$num <- row.names(m5) 
    all <- merge(m1, m2, by = "num", all = T) 
    var <- merge(all, m3, by = "num", all = T) 
    var <- merge(var, m4, by = "num", all = T)
    var <- merge(var, m5, by = "num", all = T) 
    var[is.na(var)] <- 0 
    var$num <- NULL
    var$Total_Variance <- var[,1] + var[,2] + var[,3] + var[,4] + var[,5]
    var$Proportion_Small_Effects <- var[,1]/var[,6]
    var$Proportion_Medium_Effects <- var[,2]/var[,6]
    var$Proportion_Large_Effects <- var[,3]/var[,6]
    var$Proportion_Extra_Large_Effects <- var[,4]/var[,6]
    var$Proportion_Extra_Plus_Large_Effects <- var[,5]/var[,6]
    output$Proportion_Small_Effects <- mean(var$Proportion_Small_Effects) 
    output$Proportion_Medium_Effects <- mean(var$Proportion_Medium_Effects) 
    output$Proportion_Large_Effects <- mean(var$Proportion_Large_Effects) 
    output$Proportion_Extra_Large_Effects <- mean(var$Proportion_Extra_Large_Effects)
    output$Proportion_Extra_Plus_Large_Effects <- mean(var$Proportion_Extra_Plus_Large_Effects)
  }
  
  # Export
  write.table(var, 
    file = paste0("summary/", gsub(".csv*", "", file), 
    "_varianceexplained_periteration.tsv"), 
	row.names = F, quote = F, sep = "\t") 
  write.table(output, 
    file = paste0("summary/", gsub(".csv*", "", file), 
    "_varianceexplained.tsv"), 
	row.names = F, quote = F, sep = "\t") 
}


# (3) Calculate Posterior Inclusion Probability (PIP) of CpGs as well as median/mean betas
##########################################################################################

cpgs <- names$Marker
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
anno <- anno[cpgs,c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]

for (i in loop) { 
  tmp <- fread(paste0("beta/", i))
  tmp <- as.data.frame(tmp)
  print(i)
  A <- gsub(".csv*", "", i)
  median <- apply(tmp, 2, median)
  median_LLCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.025))
  median_LCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.05))
  median_HCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.95))
  median_HHCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.975))
  names(median) <- cpgs
  names(median_LLCI) <- cpgs
  names(median_LCI) <- cpgs
  names(median_HCI) <- cpgs
  names(median_HHCI) <- cpgs
  median <- as.matrix(median)
  median_LCI <- as.matrix(median_LCI)
  median_LLCI <- as.matrix(median_LLCI)
  median_HCI <- as.matrix(median_HCI)
  median_HHCI <- as.matrix(median_HHCI)
  median_LCI <- cbind(median_LLCI, median_LCI)
  median_HCI <- cbind(median_HCI, median_HHCI)
  betas <- cbind(median,cbind(median_LCI,median_HCI))
  betas <- as.data.frame(betas)
  betas$CpG <- row.names(betas)
  betas <- betas[,c(6,1,2,3,4,5)]
  names(betas)[2:6] <- c("Median_Beta", "Beta_2.5", "Beta_5", "Beta_95", "Beta_97.5")
  write.table(betas, paste0("summary/", A, "_betas.tsv"), row.names = T, sep = "\t", quote = F)
  comp <- fread(paste0("./comp/", i))
  comp <- as.data.frame(comp)
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- cpgs
  betas1 <- cbind(anno, betas[,c(2,3)], pip[,2])
  print(A)
  write.table(betas1, paste0("summary/", A, "_medianbeta_pip.tsv"), row.names = F, sep = "\t", quote = F)

} 

# Calculate Mean Betas
for (i in loop) {

  tmp <- fread(paste0("beta/", i))
  tmp <- as.data.frame(tmp)
  print(i)
  A <- gsub(".csv*","",i)
  means <- apply(tmp, 2, mean)
  ses <- apply(tmp, 2, function(x){sd(x)/sqrt(length(x))})
  names(means) <- cpgs
  names(ses) <- cpgs
  means <- as.matrix(means)
  ses <- as.matrix(ses)  
  betas <- cbind(means,ses)
  betas <- as.data.frame(betas)
  betas$CpG <- row.names(betas)
  betas <- betas[,c(3,1,2)]
  names(betas)[2:3] <- c("Mean_Beta", "SE")
  comp <- fread(paste0("comp/", i))
  comp <- as.data.frame(comp)
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- cpgs
  betas1 <- cbind(anno, betas[,c(2,3)], pip[,2])
  write.table(betas1, paste0("summary/", A, "_meanbeta_pip.tsv"), row.names = F, quote = F, sep = "\t")

  # Sig CpGs
  for (j in c(0.95, 0.90, 0.85, 0.80)) {
    sig_cpgs <- betas1[betas1$PIP > j,]
    sig_cpgs <- sig_cpgs[order(-sig_cpgs$PIP),]
    write.table(sig_cpgs, paste0("summary/", A, "_meanbeta_pip", j,".tsv"), row.names = F, quote = F, sep = "\t")
  }

} 

# Significant CpGs sex agnostic (PIP > 0.95)
# betas_0.95 <- betas1[betas1$PIP > 0.95,] # 6 hits
# betas_0.90 <- betas1[betas1$PIP > 0.90,] # 11 hits
# betas_0.85 <- betas1[betas1$PIP > 0.85,] # 14 hits
# betas_0.80 <- betas1[betas1$PIP > 0.80,] # 20 hits

# Ola - Holy cow, multiline comment in R.
if (FALSE) {

	# Calculate Posterior Inclusion Probability of GROUPS of CpGs to avoid correlations being issue
	###########################################################################################################

	## Import data
	L <- list.files("./summary/", "_meanbeta_pip.tsv")
	pip_files <- lapply(paste0("./summary/", L), read.table, header = T, sep = "\t") 
	names <- as.character(L) 
	names <- gsub("_meanbeta_pip.*", "", names) 
	names(pip_files) <- names 
	pip_files <- Map(cbind, pip_files) 
	pip_top <- do.call(rbind, lapply(pip_files, function(x)x[x$PIP > 0.95,]))
	pip_top1 <- do.call(rbind, lapply(pip_files, function(x)x[x$PIP > 0.2,]))

	## Next, query all CpGs that are 2.5 kb < CpGs that have PIP > 0.2  
	list1 <- list()
	anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds") # Read in annotation dataframe 
	cpgs1 <- unique(pip_top1$Name) # Establish CpGs to be queried 
	for (i in 1:length(cpgs1)) { 
	# Extract the CpG of interest 
	tmp <- anno[anno$Name %in% cpgs1[i],]
	# Subset annotation to chr and basepairs +/- within 2.5kb 
	anno.tmp <- anno[anno$chr %in% tmp$chr,]
	anno.tmp <- anno.tmp[which(anno.tmp$pos <= (tmp$pos+1250)),]
	anno.tmp<- anno.tmp[which(anno.tmp$pos >= (tmp$pos-1250)),]
	anno.tmp$group <- as.character(cpgs1[i])
	# Store output 
	list1[[i]] <- anno.tmp
	# Print to denote completion 
	print(i)
	} 
	# Combine 
	list2 <- as.data.frame(do.call("rbind",list1))
	# Establish if they are the lead CpG being queried 
	list2$lead <- ifelse(list2$Name == list2$group, "lead", "other")

	# Next, ask if those CpGs have correlation coefficient >0.5 with the lead CpG in the group 
	meth <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")
	samps <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds") # Subset to samples passing QC 
	meth <- meth[,which(colnames(meth) %in% samps$Sample_Sentrix_ID)]
	# Subset to probes of interest 
	meth1 <- meth[which(row.names(meth) %in% list2$Name),] # Subset to probes of interest 
	# Set up list to store new outputs 
	keeps <- list()
	# Loop through each group to get correlation matrix 
	for(i in unique(list2$group)){ 
	# Subset to that group 
	tmp <- list2[which(list2$group %in% i),]  
	# Subset methylation data to probes 
	meth.tmp <- meth1[which(row.names(meth1) %in% tmp$Name),]
	# Correlation matrix
	cor.tmp <- cor(as.matrix(t(meth.tmp)))
	# Subset correlation matrix to lead probe 
	cor.tmp1 <- cor.tmp[which(row.names(cor.tmp)%in%tmp[tmp$lead%in%"lead","Name"]),]
	# Remove the probe being queried as it will have a correlation of 1 with itself 
	cor.tmp2 <- cor.tmp1[which(!names(cor.tmp1)%in%tmp[tmp$lead%in%"lead","Name"])]
	# Identify those (if any) with absolute correlation coefficient > 0.5 with probe being queried 
	highcor <- cor.tmp2[abs(cor.tmp2)>=0.5]
	lowcor <- cor.tmp2[abs(cor.tmp2)<0.5]
	# Separate probes into three groups - the lead CpG, those to be kept based on correlation and those to be removed
	tmp$keep <- NA
	tmp[which(tmp$Name %in% names(highcor)),"keep"] <- "keep"
	tmp[which(tmp$Name %in% names(lowcor)),"keep"] <- "remove"
	tmp[tmp$lead%in%"lead","keep"] <- "lead"
	# Store the outputs in a new file so that we can use it again for the PIP stage 
	keeps[[i]] <- tmp
	# Print to denote completion
	print(which(unique(list2$group) %in% i))
	} 

	# Combine outputs and remove those that failed the correlation stage 
	keeps1=as.data.frame(do.call("rbind",keeps))
	keeps1=keeps1[which(keeps1$keep %in% c("lead","keep")),]

	## Now, we need to get the group PIP for all of the assigned groups 
	# Loop through groups 
	for(i in unique(keeps1$group)){ 
	# Find traits which this was significant for 
	traits <- pip_top1[pip_top1$Marker %in% i,"Trait"]
	traits <- gsub("_processed", "", traits)
	# Define which files we will loop through for this CpG 
	ind <- which(names(pip_files) %in% traits)
	for(j in ind){ 
		# Get file 
		file=pip_files[[j]]
		# Subset to that group 
		tmp=keeps1[which(keeps1$group %in% i),] 
		# Calculate group PIP 
		tmp$group_PIP=sum(file[which(file$Marker %in% tmp$Name),"PIP"])
		tmp$trait=names(pip_files)[j]
		# Tidy up file 
		tmp1 <- tmp[,c("Name", "chr", "pos", "UCSC_RefGene_Name", "group", "lead", "keep", "group_PIP", "trait")]
		write.csv(tmp1, paste0("summary/grouppip_", i, "_", names(pip_files)[[j]],".csv"),row.names=F)
	}
	# Print to denote completion 
	print(which(unique(keeps1$group) %in% i))
	} 


	## Combine group PIP results 
	setwd("/Cluster_Filespace/Marioni_Group/Rob/BayesR/Group_PIP/")
	# Combine files 
	pips=list.files(".",".")
	pips1=as.data.frame(do.call("rbind", lapply(pips, read.csv, header = T)))
	# Only subset to those where group PIP > 0.8
	pips2=pips1[pips1$group_PIP >= 0.8,]
	pips3=pips2[pips2$lead%in%"lead",]
	# Tabulate results 
	table(pips3$trait)
}