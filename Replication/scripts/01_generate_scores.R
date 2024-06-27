#########################################################################
################### SCRIPT 1 : CALCULATE SCORES #1-4   ##################
#########################################################################
# Check if missing and install if so 

list.of.packages <- c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      "minfi", "htmlTable", "rmarkdown", "devtools", "EpiSmokEr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

if ("EpiSmokEr" %in% new.packages) {
  library(devtools)
  install_github("sailalithabollepalli/EpiSmokEr")
}
library(EpiSmokEr)

# Set working directory
# setwd("<path_to_replication_folder>")

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
# Of note, this example is a very truncated dataframe and most CpGs will be missing, please input your typical methylation dataframe here (all available CpGs)
data = readRDS("examples/methylation_df_example.RDS") 


## Of note, you will need CpGs as rows and IDs as columns before progressing ## 
# Check if data needs to be transposed 
if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-as.data.frame(t(data))
}

#####################################
### 3. Read in phenotype object #####
#####################################

# Read in smoking data 
# ensure that smoking in pack years is (log+1)-transformed prior to this with outliers removed, please refer to analysis plan for further info 
pheno = read.csv("examples/pheno_example.csv")

#####################################################################
### 4. Read in predictors and their weights, generate EpiSmokEr #####
#####################################################################

# Read in weights for scores 
cpgs <- read.csv("data/Predictors_long.csv", header = T) 
sites = read.csv("data/EpiSmoker_smokingScore_sites.tsv", sep='\t')
sites = c(sites$CpG, cpgs$CpG)

#################################################
### 5. Quality control prior to projections #####
#################################################

# Subset methylation dataframe to sites for prediction 
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
# coef[is.na(coef)] <- 0

# Convert NAs to mean values 
coef <- as.data.frame(t(apply(coef,1,function(x) na_to_mean(x))))


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

result <- epismoker(dataset=coef, samplesheet = pheno, method = "SSc")

# Prepare final file 
out$ID=row.names(out)
# Format 
out=out[,c(ncol(out),(1:(ncol(out)-1)))]
# merge with EpiSmoker
out = merge(out, result, by.x="ID", by.y="SampleName")
colnames(out)[5] = "EpiSmokEr"
# Merge with phenotype file
out=merge(out, pheno,by="ID",all.y=T)

# Write out to pass to next script 
write.csv(out, "results/smoking_scores.csv", row.names=F)

