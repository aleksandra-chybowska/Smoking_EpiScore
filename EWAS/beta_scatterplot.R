library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)

# Extract betas from summary statistics of all EWASs

path = '<cluster_home_dir>/TWIST/Manual_EWAS/'

source = '40_10_EWAS'
twist = read.csv(paste0(path, "pass10_40_46_full_EWAS_TWIST.csv"))

twist$study = "TWIST"
twist$chr = str_replace(twist$chr, "chr", "")

# examine non autosomal data
table(twist$chr)
mito_sex = subset(twist, chr %in% c('M', 'X', 'Y'))
table(mito_sex$p < 1e-5) 

# nothing significant, leaving only autosomes
twist$chr = as.numeric(twist$chr)
twist = subset(twist, !is.na(chr))
table(twist$chr)

ont_path = "<cluster_home_dir>/ONT/EWAS/"
ont = read.csv(paste0(ont_path, "5_40_46_full_EWAS_ONT.csv"))
ont$study = "ONT"
ont$chr = str_replace(ont$chr, "chr", "")

table(ont$chr)
mito_sex = subset(ont, chr %in% c('M', 'X', 'Y'))
table(mito_sex$p < 1e-5) 

ont$chr = as.numeric(ont$chr)
ont = subset(ont, !is.na(chr))
table(ont$chr)

# reading in EPIC data
epic = read.csv(paste0(path, 'W3_46_processed_new_annots.csv'))
epic = epic[-c(1, 10, 11, 12)]
table(epic$chromosome)
epic[epic$chromosome == "chr19_GL949749v2_alt", ] 
# not significant, fine to remove

epic$chromosome = str_replace(epic$chromosome, "chr", "")
epic$chromosome = as.numeric(epic$chromosome)
epic = subset(epic, !is.na(chromosome))

cols = c("cpg", "chr", "pos", "beta", "se", "tstat", "pval", "study")
colnames(twist) = cols
colnames(epic) = cols
colnames(ont) = cols

epic$pos = epic$pos + 1 
twist$chrpos = paste0(twist$chr, '-', twist$pos)
epic$chrpos = paste0(epic$chr, '-', epic$pos)
ont$chrpos = paste0(ont$chr, '-', ont$pos)

all = merge(twist, epic, by="chrpos", all = T)
all = merge(all, ont, by="chrpos", all = T)
saveRDS(all, "<cluster_home_dir>/ONT/EWAS/all_5_40_46_new_annots_with_beta.RDS", compress=F)

all = readRDS("<cluster_home_dir>/ONT/EWAS/all_5_40_new_annots_with_beta.RDS")
select_col = function(row, col1, col2, col3) {
  return(dplyr::first(na.omit(c(row[[col1]], row[[col2]], row[[col3]]))))
}

df = data.frame("CpG" = apply(all, 1, select_col, "cpg.y", "cpg.x", "cpg"), 
                "Chromosome" = apply(all, 1, select_col, "chr.x", "chr.y", "chr"),   
                "Position" = as.numeric(apply(all, 1, select_col, "pos.x", "pos.y", "pos")), 
                "Beta_EPIC" = all$beta.y,
                "TWIST_Beta" = all$beta.x,
                "Beta_ONT" = all$beta
)
head(df)

setwd("<cluster_home_dir>/ONT/EWAS/")
saveRDS(df, "40_10_40_5_46_plotting_ds_new_annots_betas.RDS", compress=F)

# For the CpG sites available in both platforms, 
# an additional scatter plot coupled with a regression slope for contrasting 
# the beta estimates from two platforms would be very helpful.#

# Find CpG sites present on all three platforms
# 1. If it is present on the Epic array, its CpG name will start from cg - we need only those
# 2. It can't have NA in other columns

library(gridExtra)
library(ggpointdensity)
library(tidyverse)
library(ggplot2)
library(ggpubr)          # For stat_cor to display Pearson correlation
library(patchwork)
library(cowplot)
library(data.table)

setwd("<cluster_home_dir>/ONT/EWAS/")

all = readRDS("all_5_40_46_new_annots_with_beta.RDS")
# all = readRDS("40_10_40_5_46_plotting_ds_new_annots_betas.RDS")
colnames(all) = c("chrpos", 
 "CpG_TWIST", "Chr_TWIST", "Pos_TWIST", "Beta_TWIST", "SE_TWIST", "Tstat_TWIST", "P_TWIST", "Study_TWIST",
 "CpG_EPIC", "Chr_EPIC", "Pos_EPIC", "Beta_EPIC", "SE_EPIC", "Tstat_EPIC", "P_EPIC", "Study_EPIC",  
 "CpG_ONT", "Chr_ONT", "Pos_ONT", "Beta_ONT", "SE_ONT", "Tstat_ONT", "P_ONT", "Study_ONT")

df = all[c("chrpos", "CpG_TWIST", "Beta_TWIST", "P_TWIST", 
                     "CpG_EPIC", "Beta_EPIC", "P_EPIC", 
                     "CpG_ONT", "Beta_ONT", "P_ONT")]

df_ONT_0.05 = subset(df, P_ONT <0.05)

# do we need to subset them somehow?
# I might run a filter at P<0.05 for at least one of the EWASs 
# otherwise we'll just end up with a big blob

df_all_subset = subset(df, P_ONT < 0.05 | P_TWIST < 0.05 | P_EPIC < 0.05)

dim(df)
dim(df_ONT_0.05)
head(df_ONT_0.05)

model_common_twist_ont <- lm(Beta_ONT ~ Beta_TWIST, data = df_all_subset)
r_squared_common_twist_ont <- summary(model_common_twist_ont)$r.squared #0.04216778
r_common_twist_ont <- cor(df_all_subset$Beta_ONT, df_all_subset$Beta_TWIST,
  method = "pearson", use="complete.obs")

model_common_epic_twist <- lm(Beta_EPIC ~ Beta_TWIST, data = df_all_subset)
r_squared_common_epic_twist <- summary(model_common_epic_twist)$r.squared # 0.1631512
r_common_epic_twist <- cor(df_all_subset$Beta_EPIC, df_all_subset$Beta_TWIST,
  method = "pearson", use="complete.obs")

model_common_epic_ont <- lm(Beta_EPIC ~ Beta_ONT, data = df_all_subset)
r_squared_common_epic_ont <- summary(model_common_epic_ont)$r.squared #0.1377789 
r_common_epic_ont <- cor(df_all_subset$Beta_EPIC, df_all_subset$Beta_ONT,
  method = "pearson", use="complete.obs")

colors <- c("#0c8bb4", "#44c2f2", "#ff972f", "#ffdcbf")

# Calculate the position for stat_cor label (top left corner of each plot)
label_x_pos <- -0.2  # Positions label at the far left of the plot
label_y_pos <- 0.2   # Positions label at the top of the plot
x_limit <- c(-0.2, 0.2)
y_limit <- c(-0.2, 0.2)

### With Linear regression 

# Create each scatter plot with ideal line, Pearson r, no legend, and equal coordinates
plot1 <- ggplot(df_all_subset, aes(x = Beta_TWIST, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +  # Custom colors, no legend
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", x = label_x_pos, y = label_y_pos, 
    label = bquote(r == .(sprintf("%.2f", r_common_twist_ont))),
    hjust = 0) +
  labs(x = "Beta TWIST", y = "Beta ONT") +
  xlim(x_limit) +
  ylim(y_limit) +
  coord_equal() +  # Ensures a 1:1 aspect ratio
  theme_minimal()

plot2 <- ggplot(df_all_subset, aes(x = Beta_EPIC, y = Beta_TWIST)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", x = label_x_pos, y = label_y_pos, 
    label = bquote(r == .(sprintf("%.2f", r_common_epic_twist))),
    hjust = 0) +
  labs(x = "Beta EPIC", y = "Beta TWIST") +
  xlim(x_limit) +
  ylim(y_limit) +
  coord_equal() +
  theme_minimal()

plot3 <- ggplot(df_all_subset, aes(x = Beta_EPIC, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate("text", x = label_x_pos, y = label_y_pos, 
    label = bquote(r == .(sprintf("%.2f", r_common_epic_ont))),
    hjust = 0) +
  labs(x = "Beta EPIC", y = "Beta ONT") +
  xlim(x_limit) +
  ylim(y_limit) +
  coord_equal() +
  theme_minimal()

pdf("Combined_scatterplots_with_LM_equal_coord_ABC_r_46_samps.pdf", width = 8, height = 8)
  plot_grid(plot1, plot2, plot3, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
dev.off()

#### OBSOLETE

# Create each scatter plot with ideal line, Pearson r, no legend, and equal coordinates
plot1 <- ggplot(common_twist_ont, aes(x = Beta_TWIST, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +  # Custom colors, no legend
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Ideal correlation line
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, hjust = 0, vjust = 1) +
  labs(x = "Beta TWIST", y = "Beta ONT") +
  xlim(x_limit) +
  ylim(y_limit) +
  coord_equal() +  # Ensures a 1:1 aspect ratio
  theme_minimal()

plot2 <- ggplot(common_epic_twist, aes(x = Beta_EPIC, y = Beta_TWIST)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, hjust = 0, vjust = 1) +
  labs(x = "Beta EPIC", y = "Beta TWIST") +
  xlim(x_limit) +
  ylim(y_limit) +
  coord_equal() +
  theme_minimal()

plot3 <- ggplot(common_epic_ont, aes(x = Beta_EPIC, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  stat_cor(method = "pearson", label.x = -Inf, label.y = Inf, hjust = 0, vjust = 1) +
  labs(x = "Beta EPIC", y = "Beta ONT") +
  xlim(x_limit) +
  ylim(y_limit) +
  coord_equal() +
  theme_minimal()

# Combine plots into a 1x3 grid
combined_plot <- plot1 + plot2 + plot3 + plot_layout(ncol = 3)

pdf("Combined_scatterplots_with_ideal_line_equal_coord.pdf", width = 12, height = 4)
print(combined_plot)
dev.off()

pdf("Combined_scatterplots_with_ideal_line_equal_coord_ABC.pdf", width = 8, height = 8)
plot_grid(plot1, plot2, plot3, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
dev.off()

#############################

pdf("ONT_vs_TWIST_beta_scatterplot.pdf", width = 5, height = 5)

# Plot
# Create the ggplot scatterplot with a line of best fit
ggplot(common_twist_ont, aes(x = TWIST_Beta, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors) +  # Scatterplot points in blue
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Line of best fit in red, without confidence interval shading
  labs(
    x = "Beta TWIST",
    y = "Beta ONT"
  ) +
  theme_minimal()  # Optional theme for a clean look
dev.off()

pdf("EPIC_vs_TWIST_beta_scatterplot.pdf", width = 5, height = 5)
# Plot
# Create the ggplot scatterplot with a line of best fit
ggplot(common_epic_twist, aes(x = Beta_EPIC, y = TWIST_Beta)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors) + # Scatterplot points in blue
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Line of best fit in red, without confidence interval shading
  labs(
    x = "Beta EPIC",
    y = "Beta TWIST"
  ) +
  theme_minimal()  # Optional theme for a clean look
dev.off()

pdf("EPIC_vs_Beta_ONT_scatterplot.pdf", width = 5, height = 5)
# Plot
# Create the ggplot scatterplot with a line of best fit
ggplot(common_epic_ont, aes(x = Beta_EPIC, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors) +  # Scatterplot points in blue
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Line of best fit in red, without confidence interval shading
  labs(
    x = "Beta EPIC",
    y = "Beta ONT"
  ) +
  theme_minimal()  # Optional theme for a clean look
dev.off()

# grid?

library(ggplot2)
library(ggpointdensity)  # For geom_pointdensity
library(ggpubr)          # For stat_cor to display Pearson correlation
library(patchwork)       # For arranging plots in a grid

# Define the custom color palette
colors <- c( "#0c8bb4", "#44c2f2", "#ff972f", "#ffdcbf")

# Calculate the position for stat_cor label (top left corner of each plot)
label_x_pos <- -0.3  # Positions label at the far left of the plot
label_y_pos <- 0.3   # Positions label at the top of the plot
x_limit <- c(-0.3, 0.3)
y_limit <- c(-0.3, 0.3)

# Create each scatter plot with a line of best fit, Pearson r, and without legend
plot1 <- ggplot(common_twist_ont, aes(x = Beta_TWIST, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +  # Custom colors, no legend
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(method = "pearson", label.x = label_x_pos, label.y = label_y_pos, hjust = 0, vjust = 1) +
  labs(x = "Beta TWIST", y = "Beta ONT") +
  xlim(x_limit) +
  ylim(y_limit) +
  theme_minimal()

plot2 <- ggplot(common_epic_twist, aes(x = Beta_EPIC, y = Beta_TWIST)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(method = "pearson", label.x = label_x_pos, label.y = label_y_pos, hjust = 0, vjust = 1) +
  labs(x = "Beta EPIC", y = "Beta TWIST") +
  xlim(x_limit) +
  ylim(y_limit) +
  theme_minimal()

plot3 <- ggplot(common_epic_ont, aes(x = Beta_EPIC, y = Beta_ONT)) +
  geom_pointdensity() +
  scale_color_gradientn(colors = colors, guide = "none") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  stat_cor(method = "pearson", label.x = label_x_pos, label.y = label_y_pos, hjust = 0, vjust = 1) +
  labs(x = "Beta EPIC", y = "Beta ONT") +
  xlim(x_limit) +
  ylim(y_limit) +
  theme_minimal()

# Combine plots into a 1x3 grid
combined_plot <- plot1 + plot2 + plot3 + plot_layout(ncol = 3)

# Save to a single PDF file
pdf("Combined_scatterplots.pdf", width = 15, height = 5)
print(combined_plot)
dev.off()