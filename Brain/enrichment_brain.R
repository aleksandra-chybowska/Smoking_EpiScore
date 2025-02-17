library(missMethyl)

setwd("<cluster_home_dir>/Smoking/Brain_vs_Blood/")

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
results_hc = readRDS("results_hc.RDS")
results_pfc = readRDS("results_pfc.RDS")
results_ba17 = readRDS("results_ba17.RDS")
results_ba24 = readRDS("results_ba24.RDS")
results_ba2021 = readRDS("results_ba2021.RDS")

threshold = 1e-5

significant_hc = subset(results_hc, hc_p < threshold)
significant_pfc = subset(results_pfc, pfc_p < threshold)
significant_ba17 = subset(results_ba17, ba17_p < threshold)
significant_ba24 = subset(results_ba24, ba24_p < threshold)
significant_ba2021 = subset(results_ba2021, ba2021_p < threshold)

write.csv(significant_hc, "significant_hc.csv", row.names = F)
write.csv(significant_pfc, "significant_pfc.csv", row.names = F)
write.csv(significant_ba17, "significant_ba17.csv", row.names = F)
write.csv(significant_ba24, "significant_ba24.csv", row.names = F)
write.csv(significant_ba2021, "significant_ba2021.csv", row.names = F)


nrow(significant_hc) + nrow(significant_pfc) + nrow(significant_ba17) + nrow(significant_ba24) + nrow(significant_ba2021) # 22

significant_all = c(significant_hc$hc_cg, 
                    significant_pfc$pfc_cg, 
                    significant_ba17$ba17_cg, 
                    significant_ba24$ba24_cg, 
                    significant_ba2021$ba2021_cg)
write.csv(significant_all, "significant_all.csv", row.names = F)
duplicated(significant_all) # all unique


# Randomly select 1000 CpGs to be significantly differentially methylated
sigcpgs <- significant_all
# All CpG sites tested
allcpgs <- results_hc$hc_cg

# GO testing with prior probabilities taken into account
# Plot of bias due to differing numbers of CpG sites per gene
gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO", 
              plot.bias = TRUE, prior.prob = TRUE, anno = ann)

top_gst <- topGSA(gst, n=55)
write.csv(top_gst, 'Gene_enrichment.csv')
