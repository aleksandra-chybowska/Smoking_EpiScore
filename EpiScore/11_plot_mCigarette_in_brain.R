library("foreign")
library("tidyverse")
library("gridExtra")
library("grid")
 
get_title = function(index) {
    if (index == "hc") {
        return("Hippocampus (BA35)") 
    } else if (index == "pfc") {
        return("Dorsolateral Prefrontal Cortex (BA46)")
    } else if (index == "ba17") {
        return("Primary Visual Cortex (BA17)")
    } else if (index == "ba24") {
        return("Anterior Cingulate Cortex (BA24)")
    } else if (index == "ba2021") {
        return("Ventral/Lateral Inferior Temporal Cortex (BA20/21)")
    }
}

setwd("/Cluster_Filespace/Marioni_Group/Anna/Brain_methylation/")
dat = readRDS("Normalised_Beta_Values_danet.rds")
descriptions = readRDS("descriptions.rds") # According to this file, this patient does not have HC data: LBC360666

setwd("/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore_tissues/LBC")
dat = t(dat)
identical(rownames(descriptions), rownames(dat))

dat = merge(descriptions, dat, by="row.names")

# I need to change this into pack years from my previous files for testing against this metric
ph36 = read.spss("/Cluster_Filespace/Marioni_Group/Riccardo/Cognitive_EWAS_GS_LBC/LBC1936_BloodBasedEWASofCognitiveAbilities_RM_06MAY2021.sav", to.data.frame=T, use.value.labels=F)
ph361 = ph36[,c(1,12)]
names(ph361)[1] = "LBC_ID"

d = merge(ph361, dat, by="LBC_ID")
d1 = d[-grep("_R", d$Sample_Name),] #  69 807184

tissues = list(
    hc = d1[grep("HC", d1$Sample_Name),], # smokcat_w1 missing for one person with hc
    pfc = d1[grep("BA46", d1$Sample_Name),],
    ba17 = d1[grep("BA17", d1$Sample_Name),],
    ba24 = d1[grep("BA24", d1$Sample_Name),],
    ba2021 = d1[grep("BA20", d1$Sample_Name),]
)

run = "W1_W3_W4_1256_methylation_adjusted_set_correct_target"
results = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final/"

## Predictions in LBC
####################################################################

plots = list()

for (i in names(tissues)){
    tissue = tissues[[i]]

    beta = read.csv(paste0(results, "weights_", run, ".csv"))
    rownames(beta) = beta$CpG_Site
    beta = beta[2:nrow(beta),] #1165    3

    # Match CpGs    
    beta = beta[which(rownames(beta) %in% colnames(tissue)),]
    brain_dnam = NA
    brain_dnam = tissue[,rownames(beta)]

    # Scale
    brain_dnam = scale(brain_dnam)

    # Replace missing values with 0
    brain_dnam[is.na(brain_dnam)] = 0

    pred_36 = brain_dnam %*% beta$Coefficient
    pred_36 = as.data.frame(pred_36)
    names(pred_36) = c("py_pred")
    pred_36["smokcat_w1"] = factor(tissue$smokcat_w1, labels=c("Never", "Ex", "Current"))

    plots[[i]] = ggplot(pred_36, aes(x=as.factor(smokcat_w1), y=py_pred, color=as.factor(smokcat_w1))) + # fill=name allow to automatically dedicate a color for each group
        scale_color_manual(values=c("#C94905", "#218F60", "#155DA2")) +
        geom_point(size=0.7, alpha=0.9,show.legend = FALSE, position=position_jitter(w = 0.15, h = 0)) +
        ggtitle(get_title(i)) +
        theme(plot.title = element_text(size=10),
        axis.title = element_text(size=10)) +
        xlab("") +
        ylab("") +
        ylim(c(-1.5, 1.5)) 
}

y.grob <- textGrob("mCigarette", gp=gpar(fontsize=12), rot=90)

ggsave(file = "Grid_mCigarette_in_brain.pdf", arrangeGrob(grobs = plots, ncol = 2, left = y.grob))

pdf(file="Grid_mCigarette_in_brain.pdf", width = 4, height = 4)
plots$hc
# grid.arrange(arrangeGrob(p3, p4, top = grid::textGrob("Blood", x = 0.5, hjust = 0)),
#              arrangeGrob(p2, p1, top = grid::textGrob("Hippocampus", x = 0.35, hjust = 0)), ncol = 2)
dev.off()
