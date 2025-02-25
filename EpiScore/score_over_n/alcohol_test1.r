#########################################################################
################### SCRIPT 1 : CALCULATE SCORES #1-4   ##################
#########################################################################
# Check if missing and install if so 
#merged = readRDS("lbc_comp_merged.RDS")

# list.of.packages <- c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                       "minfi", "htmlTable", "rmarkdown", "devtools", "EpiSmokEr")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

# # Load packages
# library(IlluminaHumanMethylation450kmanifest)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(minfi)
# library(htmlTable)
# library(rmarkdown)

# if ("EpiSmokEr" %in% new.packages) {
#   library(devtools)
#   install_github("sailalithabollepalli/EpiSmokEr")
# }
# library(EpiSmokEr)

# Set working directory
wd = "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/score_over_n_alcohol/"
setwd(wd)

#####################################################
######## 1. Create functions for tidying data #######
#####################################################

# M2beta 
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

# NA to mean (mean imputation)
na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
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
rownames(target) <- target$Basename
target <- target[target["cohort"] == "LBC36" ,] # 2797, waves 1, 2, 3, and 4
target <- subset(target, WAVE==1)
data <- data[, rownames(target)] # 895 1349

# do tej pory dobrze

df = data[1:10, ]
df = as.data.frame(df)
rownames(df) = paste0("ID", 1:nrow(df))
saveRDS(df, "methylation_df_example.rds")


## Of note, you will need CpGs as rows and IDs as columns before progressing ## 
# this is the opposite form of what is expected in script_02 # 
# however, script should handle this for you # 
# Check if Data needs to be Transposed 
if(ncol(data) > nrow(data)){ # REMOVE TRUE FROM HERE
  message("It seems that individuals are rows - data will be transposed!")
  data<-as.data.frame(t(data))
}

#####################################
### 3. Read in phenotype object #####
#####################################

# Read in CRP # ensure that CRP is log-transformed prior to this with outliers removed, please refer to analysis plan for further info 
pheno = read.table("<filespace_marioni_group_dir>/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T,
 row.names = 1)
pheno = subset(pheno, cohort == "LBC36" & WAVE == 1)
pheno = pheno[c("age", "sex", "alcunitsupw", "alcohol_cat", "drinker")]
pheno$Basename = rownames(pheno)

pheno_example = pheno[1:10, ]
pheno_example$Basename = paste0("ID", 1:10)
rownames(pheno_example) = pheno_example$Basename
pheno_example = pheno_example[1:5]
saveRDS(pheno_example, "alcohol_pheno_example.rds")


#####################################################################
### 4. Read in predictors and their weights, generate EpiSmokEr #####
#####################################################################

cpgs <- NULL
files = list.files(path = wd, pattern="weights_")

for (file in files) {
  cp <- read.csv(file, header = T)
  name <- gsub("[^0-9-]", "", file)
  cp$Predictor <- as.numeric(name)
  cpgs <- rbind(cpgs, cp)
}

colnames(cpgs) = c("CpG", "Beta", "Predictor")

# EpiSmoker - separate two scripts

#################################################
### 5. Quality control prior to projections #####
#################################################

# Subset methylation dataframe to sites for prediction 
sites = cpgs$CpG
coef = data[intersect(rownames(data), sites),]

# Convert M values to Beta values 
data<-if((range(data,na.rm=T)> 1)[[2]] == "TRUE") { 
  message("Suspect that M Values are present. Converting to Beta Values")
  data = m_to_beta(data) 
} else { 
  message("Suspect that Beta Values are present");
  data
}
cols = colnames(coef)
coef = as.data.frame(t(apply(coef, 1, scale)))
colnames(coef) = cols
#coef[is.na(coef)] <- 0

# Convert NAs to mean values 
coef <- as.data.frame(t(apply(coef,1,function(x) na_to_mean(x)))) # this may cause issues here?


#################################################
### 6. Loop stage for projecting scores  ########
#################################################

# Loop through each predictor 
loop = unique(cpgs$Predictor)
# Set up output dataframe 
out <- data.frame()
# Initiate loop 
for(i in loop){ 
  # Find CpGs for predictor i 
  tmp=coef[intersect(row.names(coef), cpgs[cpgs$Predictor %in% i,"CpG"]),]
  
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  # If there is more than one CpG present for the predictor (hopefully will be), then align CpGs across files
  if(nrow(tmp_coef) > 1) { 
    # Align CpGs in your methylation dataframe with those in the predictor file 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG),]
    # Multiply methylation beta values by weight for each CpG - and take sum for each individual 
    out[colnames(coef),paste0(i, "CpGs")]=colSums(tmp_coef$Beta*tmp)
  } else {
    # If there is just one CpG, then it will already be aligned with dataframe, so can just calculate and take score 
    tmp2 = as.matrix(tmp)*tmp_coef$Beta 
    out[colnames(coef),paste0(i, "CpGs")] = colSums(tmp2)
  }
} 

#############################################################
### 7. Prepare final file and merge with pheno file  ########
#############################################################

# Prepare final file 
out$Basename=row.names(out)
# Format 
out=out[,c(ncol(out),(1:(ncol(out)-1)))]

# Merge with phenotype file
out=merge(out, pheno,by="Basename",all.y=T)

# Write out to pass to next script 
write.csv(out, "alcohol_scores_by_n.csv", row.names=F)
