library(data.table)
smoking = fread('<cluster_home_dir>/Smoking/LDSC/Summary_statistics/EUR_without_UKB/EUR_stratified_without_UKB/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz')
# dim(smoking)
# [1] 13924395       10
# min N: 1839
# median: 180448
# max N: 183196


smoking_with_UKB=fread('<cluster_home_dir>/Smoking/LDSC/Summary_statistics/EUR_stratified/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.gz')
# dim(smoking_with_UKB)
# [1] 13763312       10
# min N : 3455
# median N: 323872
# max N: 326497