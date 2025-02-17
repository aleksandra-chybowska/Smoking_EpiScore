#########################################################################
################### SCRIPT 1 : CALCULATE ALCOHOL SCORE   ################
#########################################################################

# Set working directory
setwd("<cluster_home_dir>/Smoking/Replication_Elena/test")

#####################################################
######## 1. Create functions for tidying data #######
#####################################################

# M2beta 
m_to_beta = function (val) 
{
  beta = 2^val/(2^val + 1)
  return(beta)
}

# NA to mean (mean imputation)
na_to_mean = function(methyl){
  methyl[is.na(methyl)]=mean(methyl,na.rm=T)
  return(methyl)
}

#####################################
### 2. Read in methylation object ###
#####################################

# Methylation dataframe 
data = readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds") # Of note, this example is a very truncated dataframe and most CpGs will be missing, please input your typical methylation dataframe here (all available CpGs)
data = as.data.frame(data, check.names=F)
target = read.table("<filespace_marioni_group_dir>/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T,
 row.names = 1)
target$Basename = rownames(target)
target = target[target["cohort"] == "LBC36" ,] # 2797, waves 1, 2, 3, and 4
target = subset(target, WAVE==1)
data = data[, rownames(target)] # 895 1349

#####################################
### 3. Read in phenotype object #####
#####################################

# Read in CRP # ensure that CRP is log-transformed prior to this with outliers removed, please refer to analysis plan for further info 
pheno = read.table("<filespace_marioni_group_dir>/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T,
 row.names = 1)
pheno = subset(pheno, cohort == "LBC36" & WAVE == 1)
pheno = pheno[c("age", "sex", "alcunitsupw")]
pheno$Basename = rownames(pheno)
pheno$log_alcohol = log(pheno$alcunitsupw + 1)
pdf("alc_hists.png")
  hist(pheno$log_alcohol)
  hist(pheno$alcunitsupw)
dev.off()

pheno = pheno[c("Basename", "age", "sex", "log_alcohol")]

#####################################################################
### 4. Read in predictors and their weights, generate EpiSmokEr #####
#####################################################################

# Read in weights for scores 
cpgs = read.csv("weights.csv")

#################################################
### 5. Quality control prior to projections #####
#################################################

# Subset methylation dataframe to sites for prediction 

sites = cpgs$CpG
coef = data[intersect(rownames(data), sites),]

# Convert M values to Beta values 
data=if((range(data,na.rm=T)> 1)[[2]] == "TRUE") { 
  message("Suspect that M Values are present. Converting to Beta Values")
  data = m_to_beta(data) 
} else { 
  message("Suspect that Beta Values are present");
  data
}
cols = colnames(coef)
coef = as.data.frame(t(apply(coef, 1, scale)))
colnames(coef) = cols

# Convert NAs to mean values 
coef = as.data.frame(t(apply(coef,1,function(x) na_to_mean(x))))

#################################################
### 6. Loop stage for projecting scores  ########
#################################################

# Loop through each predictor 
loop = unique(cpgs$Predictor)
# Set up output dataframe 
out = data.frame()
# Initiate loop 
for(i in loop){ 
  # Find CpGs for predictor i 
  tmp = coef[intersect(row.names(coef), cpgs[cpgs$Predictor %in% i,"CpG"]),]
  
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  # If there is more than one CpG present for the predictor (hopefully will be), then align CpGs across files
  if(nrow(tmp_coef) > 1) { 
    # Align CpGs in your methylation dataframe with those in the predictor file 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG),]
    # Multiply methylation beta values by weight for each CpG - and take sum for each individual 
    out[colnames(coef),i]=colSums(tmp_coef$Beta*tmp)
  } else {
    # If there is just one CpG, then it will already be aligned with dataframe, so can just calculate and take score 
    tmp2 = as.matrix(tmp)*tmp_coef$Beta 
    out[colnames(coef),i] = colSums(tmp2)
  }
} 

#############################################################
### 7. Prepare final file and merge with pheno file  ########
#############################################################

# Prepare final file 
out$Basename = row.names(out)
# Format 
out = out[,c(ncol(out),(1:(ncol(out)-1)))]

# Merge with phenotype file
out = merge(out, pheno,by="Basename",all.y=T)

# Write out to pass to next script 
write.csv(out, "alcohol_scores.csv", row.names=F)

################################
#### 8. Test performance  ######
################################

# Read in previous scores
pheno = read.csv("alcohol_scores.csv")
score = "Alcohol.Consumption.Score"

r = cor(pheno$log_alcohol, pheno[score], use="pairwise.complete.obs")

# Incremental DNAm R2
null = summary(lm(pheno$log_alcohol ~ pheno$age + pheno$sex))$r.squared
full = summary(lm(pheno$log_alcohol ~ pheno$age + pheno$sex + pheno[,score]))$r.squared 
w1 = round(100*(full - null), 3) 

out = data.frame("score" = score,  
                 "null_R2" = null,
                 "full_R2" = full,
                 "incremental_R2" = w1,
                 "r" = as.numeric(r))

write.csv(out, "results.csv")
