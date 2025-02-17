library(data.table)
library(dplyr)
library(jsonlite)
setwd('<datastore_home_dir>/Lab/Smoking/GWAS/')

epismoker_included = read.csv("../EpiSmoker_CPG.csv")

epismoker = 'epismoker_overlap/outfile.fastGWA'
epismoker_res = read.table(epismoker, header = T)
epismoker_res$STUDY = "Epismoker"
min_smoker_res = epismoker_res[c("SNP", "BETA", "CHR", "POS", "P", "STUDY")]
data = subset(min_smoker_res, P<1e-5)
anno = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno, "V8", "V2")
data = merge(data, anno, by.x="SNP", by.y="V8", all.x = T)

significant = subset(min_smoker_res, P<5e-8)
significant = merge(significant, anno, by.x="SNP", by.y="V8", all.x = T)
write.csv(data, 'epismoker_overlap/suggestive_rsids.csv')
write.csv(significant, 'epismoker_overlap/significant_rsids.csv')
rsids = significant$V2
suggestive_rsids = data$V2

json = toJSON(rsids)
suggestive_json = toJSON(suggestive_rsids)

writeLines(json, "epismoker_overlap/significant.json")
writeLines(suggestive_json, "epismoker_overlap/suggestive_json.json")
# curl -i -H "Content-Type: application/json" -X POST -d @significant.json http://api.godmc.org.uk/v0.1/query

significant_response = fromJSON("epismoker_overlap/significant_response.json")
overlap = intersect(epismoker_included$CpG, significant_response$cpg)

