setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/")

loop = list.files("Sigma/Methylation/", pattern = ".csv")
file = "pack_years.csv"

output = matrix(nrow = 1, ncol = 1) 
output <- as.data.frame(output) 
names(output)[1] <- "Trait" 

sigma <- read.csv(paste("Sigma/Methylation/", file, sep=""))

output$Epigenetic_Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T)
output$Epigenetic_Low_CI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
output$Epigenetic_High_CI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)


library(data.table)
betas <- fread(paste("Beta/Methylation/", file, sep=""))
betas = as.data.frame(betas)
comp <- fread(paste("Comp/Methylation/", file, sep=""))  
comp = as.data.frame(comp)

nm = gsub(".csv", "", file)
names=read.csv(paste0("CpGs/", nm, "_cpg_list.csv"),header=T)
names(comp) <- names$Marker
names(betas) <- names$Marker

list = apply(comp,1,function(x)which(!x%in% c(1,2,3)))  
  
  x <- as.matrix(0, ncol = 1, nrow = 1000) 
  x <-as.data.frame(x) 
  for(i in 1:1000){ 
    x[[i]]<-length(list[[i]]) == ncol(comp) ##### NUMBER OF PROBES/MARKERS STUDIED HERE ###### 
  } 
  
  if(length(which(x[1,] %in% "TRUE")) > 0){ 
    comp = comp[-which(x %in% "TRUE"),] 
  } else { 
    comp = comp
  } 
  
  t <- vector() 
  list = apply(comp,1,function(x)which(x %in% 1)) 
  for(i in 1:1000){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true = which(t %in% "TRUE")
  ind_false = which(t %in% "FALSE")
  list_true = list[ind_true]
  list_false = list[ind_false] 
  n = length(list_true) 
  m1_1_true <- matrix(0, ncol = 1, nrow = n)
  m1_1_true <- as.data.frame(m1_1_true) 
  m1_1_true$ind <- ind_true
  x<-vector()
  for(j in m1_1_true$ind){ 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m1= as.data.frame(x) 
  m1$x[is.na(m1$x)] <- 0 
  names(m1) <- "Variance_Small_Effects" 
  
  
  t <- vector() 
  list = apply(comp,1,function(x)which(x %in% 2)) 
  for(i in 1:1000){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true = which(t %in% "TRUE")
  ind_false = which(t %in% "FALSE")
  list_true = list[ind_true]
  list_false = list[ind_false] 
  n = length(list_true) 
  m2_true <- matrix(0, ncol = 1, nrow = n)
  m2_true <- as.data.frame(m2_true) 
  m2_true$ind <- ind_true
  x<-vector()
  for(j in m2_true$ind){ 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m2= as.data.frame(x) 
  m2$x[is.na(m2$x)] <- 0 
  names(m2) <- "Variance_Medium_Effects"
  
  t <- vector() 
  list = apply(comp,1,function(x)which(x %in% 3)) 
  for(i in 1:1000){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true = which(t %in% "TRUE")
  ind_false = which(t %in% "FALSE")
  list_true = list[ind_true]
  list_false = list[ind_false] 
  n = length(list_true) 
  m3_true <- matrix(0, ncol = 1, nrow = n)
  m3_true <- as.data.frame(m3_true) 
  m3_true$ind <- ind_true
  x<-vector()
  for(j in m3_true$ind){ 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m3= as.data.frame(x) 
  m3$x[is.na(m3$x)] <- 0 
  names(m3) <- "Variance_Large_Effects"
  
  m1$num <- row.names(m1) 
  m2$num <- row.names(m2) 
  m3$num <- row.names(m3) 
  all = merge(m1, m2, by = "num", all = T) 
  var = merge(all, m3, by = "num", all = T) 
  var[is.na(var)] <- 0 
  var$num <- NULL
  var$Total_Variance <- var[,1] + var[,2] + var[,3]
  var$Proportion_Small_Effects <- var[,1]/var[,4]
  var$Proportion_Medium_Effects <- var[,2]/var[,4]
  var$Proportion_Large_Effects <- var[,3]/var[,4]
  output$Proportion_Small_Effects <- mean(var$Proportion_Small_Effects) 
  output$Proportion_Medium_Effects <- mean(var$Proportion_Medium_Effects) 
  output$Proportion_Large_Effects <- mean(var$Proportion_Large_Effects) 

  A=gsub(".csv","",file)
  output$Trait <- A
  
  write.csv(output, file = paste("Summary/", A, "_output.csv", sep = ""), row.names = F) 

#   Trait Epigenetic_Mean_Variance_Explained Epigenetic_Low_CI Epigenetic_High_CI
# 1    NA                          0.4212081         0.3509789          0.5085692


