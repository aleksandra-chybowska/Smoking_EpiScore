library(miamiplot)
library(stringr)
library(dplyr)
library(data.table)
library(CMplot)
#View(gwas_results)
setwd('U:/Datastore/IGMM/marioni-lab/Ola/Lab/Smoking/GWAS/')
tpy = 'transformed_pack_years_overlap/outfile.fastGWA'
largest_py = 'largest_gwas_pack_years/ErzurumluogluAM_30617275_discovery-stage_meta-analysis_PackYears.txt'
epismoker = 'epismoker_overlap/outfile.fastGWA'
mode = "largest" # tpy for our results, largest for 3rd party

epismoker_res = read.table(epismoker, header = T)
epismoker_res$STUDY = "Epismoker"
min_smoker_res = epismoker_res[c("SNP", "BETA", "CHR", "POS", "P", "STUDY")]

if(mode == "tpy") {
  res = read.table(tpy, header = T)
  res$STUDY = "Transformed_pack_years"
  min_ds = res[c("SNP", "CHR", "BETA", "POS", "P", "STUDY")]
} else {
  res = read.table(largest_py, header = T)
  res$STUDY = "Pack_years_GWAS"
  min_ds = res[c("rsID", "CHROM", "EFFECT_SIZE", "POS", "PVALUE", "STUDY")]
  colnames(min_ds) = c("SNP", "CHR", "BETA", "POS", "P", "STUDY")
}

data = rbind(min_ds, min_smoker_res)
anno = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno, "V8", "V2")
data = merge(data, anno, by.x="SNP", by.y="V8", all.x = T)
data$SNP = ifelse(is.na(data$V2), data$SNP, data$V2)

rm(anno)
rm(epismoker_res)
rm(res)
gc()

data = subset(data, P < 0.05)
data = subset(data, CHR != "X")

data$POS = as.numeric(data$POS)
data$CHR = as.numeric(data$CHR)
mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(22)
mycolors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]                                                                 

plot_data = prep_miami_data(data = data, split_by = "STUDY", split_at = "Epismoker", p = "P", chr = "CHR", pos = "POS")

studyA_labels <- plot_data$upper %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0(SNP)) %>%
  select(rel_pos, logged_p, label) %>%
  arrange(desc(logged_p)) %>%
  slice(1:5)

# Study B
studyB_labels <- plot_data$lower %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0(SNP)) %>%
  select(rel_pos, logged_p, label) %>%
  arrange(desc(logged_p)) %>%
  slice(1:5)

ggmiami(data = data, 
        split_by = "STUDY", split_at = "Epismoker", p = "P", chr = "CHR", pos = "POS",
        upper_ylab = "EpiSmokEr",
        lower_ylab = "Largest GWAS",
        upper_labels_df = studyA_labels, lower_labels_df = studyB_labels,
        suggestive_line_color = "blue",
        genome_line = 5e-8,
        suggestive_line = 1e-5, 
        # upper_highlight = studyA_highlight_cpg$cpg,
        #  upper_highlight_col = "cpg"
)
