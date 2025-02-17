# original source <filespace_marioni_group_dir>/Daniel_TB_Lausanne/GS/final_smoking_GS_nr.history

load("<filespace_marioni_group_dir>/Daniel_TB_Lausanne/GS/GS_variables_for_DanielT_26Feb2018_v2.RData")
load("<filespace_marioni_group_dir>/Daniel_TB_Lausanne/GS/SNPMX_input.RData")
X<-merge(SNPMX,data[,c("Methylation_ID","pack_years","GWAS_ID")],by.x="row.names",by.y="Methylation_ID")
ids <- read.table("<filespace_marioni_group_dir>/Daniel_TB_Lausanne/GS/GWAS_plink_files/GS20K_QCpsychprotocol_SPH_04112015_rel0.025_meth.grm.id", header=F)
X <- X[X$GWAS_ID %in% ids$V2,]
rownames(X)<-X$Row.names
Y<-X$pack_years.y
names(Y) <- rownames(X)
X$pack_years.y<-NULL
X$pack_years.x<-NULL
X$Row.names<-NULL
X$GWAS_ID<-NULL
X <- X[!is.na(Y),]
Y<-Y[!is.na(Y)]
Y<-scale(log(Y+1))
X<-scale(X)
X<-as.matrix(X)
y<-as.matrix(Y)
##prepare sampler
cva<-matrix(c(0.0001,0.001,0.01,0.01,0.1,1.0),nrow=2,ncol=3,byrow=T )
G=2
groupAssign=c(rep(1,380506),rep(0,941633-380506))
length(groupAssign)
library(BayesRRcpp,lib.loc="~/Rlib/Bayes")
smk_names<-colnames(X)
save(list="smk_names",file="./output/SNxM/smk_col_names_nr.RData")
savehistory("bayesrr_smk_nr_v2.history")
set.seed(234234)
any(is.na(X))
any(is.na(Y))
X[1:10,1:10]
Y[1:10]
var(Y)
var(X[,1])
X[,"Female"]

BayesRSamplerV2Groups("./output/SNxM/C1smk_GS_nr.csv",2,20000,10000,10,X,Y,0.001,0.001,0.001,0.001,0.001,cva,2.0,groupAssign)

library(data.table)
C1<-fread("./output/SNxM/C1smk_GS_nr.csv")
mean(C1$sigmaE)
mean(C1$sigmaG[1])
mean(C1$sigmaG[2])
head(C1$sigmaG[1])
C1<-as.matrix(C1)
head(C1[,"sigmaE"])
colMeans(C1[,grep("sigma",colnames(C1)])
colMeans(C1[,grep("sigma",colnames(C1))])
tmp<-C1[,"sigmaG[2]"]/rowSums(C1[,grep("sigma",colnames(C1))])
quantile(tmp)
savehistory("final_smoking_GS_nr.history")
