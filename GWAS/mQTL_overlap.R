library(data.table)
mCigarette = fread("<datastore_home_dir>/Lab/Smoking/mCigarette_CPG.csv")$CpG
epismoker = fread("<datastore_home_dir>/Lab/Smoking/EpiSmoker_CPG.csv")$CpG
mqtls = fread("<datastore_home_dir>/Lab/Smoking/Grimage_mQTL_CPG.csv")$CpG
epismoker_vs_mqtls = intersect(epismoker,mqtls)
mCigarette_vs_mqtls = intersect(mCigarette, mqtls)
