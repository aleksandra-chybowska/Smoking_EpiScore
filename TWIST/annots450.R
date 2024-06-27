if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annots = as.data.frame(annotation.table)
cpgs = c("cg09558034", "cg26381592", "cg10372485", "cg10892497", "cg08545123", 
         "cg23361916", "cg04763340", "cg07040013", "cg21292796", "cg08258765", 
         "cg00956639", "cg05265143", "cg27337459", "cg14500945", "cg13708832", 
         "cg07980535", "cg00173309", "cg17087373", "cg03425918", "cg15364032", 
         "cg25345178", "cg10395330")
sites = subset(annots, Name %in% cpgs)
min_sites = sites[c(1:4, 19, 24, 26, 32)]
min_sites = min_sites[match(cpgs, min_sites$Name),]
View(min_sites)
write.csv(min_sites, "U:/Datastore/IGMM/marioni-lab/Ola/Lab/Smoking/Brain/annots.csv")
saveRDS(min_sites, "U:/Datastore/IGMM/marioni-lab/Ola/Lab/Smoking/Brain/annots.RDS")
