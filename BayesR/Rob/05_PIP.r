library("data.table")
library("foreach")
library("doParallel")

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/PIP/Methylation/")
loop = list.files("../../Comp/Methylation/", pattern = ".csv")
library(data.table) 

i = "pack_years.csv"
nm = gsub(".csv", "", i)
names=read.csv(paste0("../../CpGs/", nm, "_cpg_list.csv"),header=T)
comp <- fread(paste("../../Comp/Methylation/", i, sep=""))  
comp<-as.data.frame(comp) 
pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
pip <- as.data.frame(reshape2::melt(unlist(pip)))
pip <- setDT(pip, keep.rownames = TRUE) 

names(pip) <- c("Marker", "PIP") 
  pip$Marker <- names$Marker
  pip$Trait <- nm 
  # Save out file 
write.csv(pip, file = paste0(nm, "_pip.csv", sep = ""), row.names = F) 
subset(pip, PIP > 0.95)
#        Marker   PIP      Trait
#  1: cg14371731 1.000 pack_years
#  2: cg00475490 1.000 pack_years
#  3: cg14391737 1.000 pack_years
#  4: cg22107533 0.991 pack_years
#  5: cg21163477 1.000 pack_years
#  6: cg21911711 1.000 pack_years
#  7: cg15159987 1.000 pack_years
#  8: cg02517189 0.995 pack_years
#  9: cg06644428 1.000 pack_years
# 10: cg23079012 1.000 pack_years
# 11: cg23288337 1.000 pack_years
# 12: cg25530474 1.000 pack_years
# 13: cg18754985 1.000 pack_years
# 14: cg23480021 1.000 pack_years
# 15: cg19859270 0.998 pack_years
# 16: cg06760111 0.992 pack_years
# 17: cg05575921 1.000 pack_years
# 18: cg04135110 1.000 pack_years
# 19: cg00731338 1.000 pack_years
# 20: cg25648203 1.000 pack_years
# 21: cg26703534 0.991 pack_years
# 22: cg14466441 0.995 pack_years
# 23: cg11207515 1.000 pack_years
# 24: cg14335029 0.989 pack_years
#         Marker   PIP      Trait


# # > nm
# # [1] "pack_years"

