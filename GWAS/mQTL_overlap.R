library(data.table)
mCigarette = fread("U:/Ola/Lab/Smoking/mCigarette_CPG.csv")$CpG
epismoker = fread("U:/Ola/Lab/Smoking/EpiSmoker_CPG.csv")$CpG
mqtls = fread("U:/Ola/Lab/Smoking/Grimage_mQTL_CPG.csv")$CpG
epismoker_vs_mqtls = intersect(epismoker,mqtls)
mCigarette_vs_mqtls = intersect(mCigarette, mqtls)
