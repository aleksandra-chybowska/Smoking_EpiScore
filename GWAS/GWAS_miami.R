library(miamiplot)
library(stringr)
library(dplyr)
library(data.table)

#View(gwas_results)
setwd('/Users/shirin/Documents/Edinburgh/Lab/Smoking/GWAS/')
tpy = 'transformed_pack_years_overlap/outfile.fastGWA'
res = read.table(tpy, header = T)
epismoker = 'epismoker_overlap/outfile.fastGWA'
epismoker_res = read.table(epismoker, header = T)

res$STUDY = "Transformed_pack_years"
epismoker_res$STUDY = "Epismoker"

min_ds_tpy = res[c("SNP", "CHR", "BETA", "POS", "P", "STUDY")]
min_smoker_res = epismoker_res[c("SNP", "BETA", "CHR", "POS", "P", "STUDY")]

data = rbind(min_ds_tpy, min_smoker_res)
anno = fread('anno_full.csv', stringsAsFactors=F, header=F)
anno = select(anno, "V8", "V2")
data = merge(data, anno, by.x="SNP", by.y="V8")
data$SNP = data$V2

rm(anno)
rm(epismoker_res)
rm(res)
gc()

data = subset(data, P < 0.05)
data$POS = as.numeric(data$POS)

mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(22)
mycolors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]                                                                 

plot_data = prep_miami_data(data = data, split_by = "STUDY", split_at = "Epismoker", p = "P", chr = "CHR", pos = "POS")

studyA_labels <- plot_data$upper %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0(SNP, "\n", BETA)) %>%
  select(rel_pos, logged_p, label) %>%
  arrange(desc(logged_p)) %>%
  slice(1:5)

# Study B
studyB_labels <- plot_data$lower %>%
  group_by(chr) %>%
  arrange(desc(logged_p)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0(SNP, "\n", BETA)) %>%
  select(rel_pos, logged_p, label) %>%
  arrange(desc(logged_p)) %>%
  slice(1:5)

ggmiami(data = data, 
        split_by = "STUDY", split_at = "Epismoker", p = "P", chr = "CHR", pos = "POS",
        upper_ylab = "EpiSmokEr",
        lower_ylab = "Transformed pack years",
        upper_labels_df = studyA_labels, lower_labels_df = studyB_labels,
        suggestive_line_color = "blue"
        # upper_highlight = studyA_highlight_cpg$cpg,
        #  upper_highlight_col = "cpg"
)
