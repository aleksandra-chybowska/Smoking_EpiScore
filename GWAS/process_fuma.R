library(tidyverse)

setwd('<datastore_home_dir>/Lab/Smoking/GWAS/FUMA_grimage')

annots = fread('<local_home_dir>/Desktop/anno_full.csv')
annots = annots[,c("V2", "V8")]
colnames(annots) = c("rsID", "GSID")

gwas_catalog = read.csv('gwascatalog.txt', sep='\t')
lead = read.csv('leadSNPs.txt', sep='\t')
lead = rename(lead, GSID = rsID)
lead = merge(lead, annots, all.x=T, by="GSID")
gwas_catalog = subset(gwas_catalog, SNPs %in% lead$rsID)
length(unique(gwas_catalog$SNPs))

write.csv(gwas_catalog, "gwascatalog_lead_loci.csv", row.names = F)

# rsid: [ids]
json = toJSON(lead$rsID)
writeLines(json, "../goDMC/lead_ids.json")
# curl -i -H "Content-Type: application/json" -X POST -d @lead_ids.json http://api.godmc.org.uk/v0.1/query > results.json

result = fromJSON("../goDMC/results.json")
result = as.data.frame(result)

write.csv(result, "../goDMC/results.csv", row.names = F)

godmc = read.csv( "../goDMC/results.csv")


