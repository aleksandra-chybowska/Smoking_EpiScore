#library(devtools)
#install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
#BiocManager::install("minfi")

library("tidyverse")
library("minfi")
library("IlluminaHumanMethylationEPICanno.ilm10b5.hg38")

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
saveRDS(annEPIC, "/Cluster_Filespace/Marioni_Group/Ola/TWIST/Annotations/IlluminaHumanMethylationEPICanno.ilm10b5.hg38.RDS", compress = F)

otherEPIC <- read.csv("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Annotations/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
ourEPICpath <-"/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/EPIC.hg38.manifest.tsv"
ourEPIC <- read.table(ourEPICpath, header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")

# according to TWIST bed cg14817997 maps to  chr1:10524-10526
subset(annEPIC, Name == "cg14817997") # chr1      10524   10526
subset(otherEPIC, Name == "cg14817997") # chr1    10524   10526 
subset(ourEPIC, Probe_ID == "cg14817997") # chr1  181054  181056 <= wrong

# AHRR
subset(annEPIC, Name == "cg05575921") # chr5      373262      373264
subset(otherEPIC, Name == "cg05575921") # chr5     373262   373264 
subset(ourEPIC, Probe_ID == "cg05575921") # chr5  373262  373264 <= OK

# how many mismatches?
otherEPICmin = otherEPIC[c("Name", "Start_hg38", "End_hg38")]
ourEPICmin = ourEPIC[c("Probe_ID", "CpG_beg", "CpG_end")]
df = merge(otherEPICmin, ourEPICmin, by.y="Probe_ID", by.x = "Name") # 865918
wrong = subset(df, Start_hg38 != CpG_beg) #  1436
# OK, will use otherEPIC from now on!

# lets find out how many sites from EPIC are available on TWIST & ONT
# TWIST
library(data.table)

setwd("/Cluster_Filespace/Marioni_Group/Ola/TWIST/Annotations/")
bed = fread("covered_targets_Twist_Methylome_hg38_annotated_separated.bed")
colnames(bed) = c("chr", "start", "stop", "cpg")

epic = read.csv("infinium-methylationepic-v-1-0-b5-manifest-file-from-Illumina.csv")
probes = epic$Name
overlap = intersect(bed$cpg, probes) # 821526/866554 = 0.9480379
saveRDS(overlap, "overlap_between_twist_and_epic.RDS")

hits = readRDS("../Manual_EWAS/twist/40_10_circos_data.RDS")
df = subset(hits, TWIST < 1e-5)
head(df)
on_illumina = subset(df, CpG %in% overlap) # 10
on_illumina_2 = subset(df, CpG %like% "cg") # 10

######################################
## new file for ONT
setwd("/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS/")
hits_ont = readRDS("40_10_40_5_plotting_ds_new_annots_46_samples.RDS")
df_ont = subset(hits_ont, TWIST < 1e-5)

on_illumina_twist_2 = subset(df_ont, CpG %in% overlap) #10

df_ont = subset(hits_ont, ONT < 1e-5)
on_illumina_ont = subset(df_ont, CpG %like% "cg") #3

######################################
## overlap

setwd("/Cluster_Filespace/Marioni_Group/Ola/ONT/EWAS")
ewas = readRDS("40_10_40_5_plotting_ds_new_annots_46_samples.RDS")
t_top = 3.6e-8
t = 1e-5

top_EPIC = subset(ewas, EPIC < t_top) # 15
top_TWIST = subset(ewas, TWIST < t_top) # 33
top_ONT = subset(ewas, ONT < t_top) # 9

suggestive_EPIC = subset(ewas, EPIC < t) # 42
suggestive_TWIST = subset(ewas, TWIST < t) # 102
suggestive_ONT = subset(ewas, ONT < t) # 63


# szukamy tego co znalazl EPIC, a czego nie znalazly inne technologie
overlap = intersect(suggestive_EPIC$CpG, probes) # 821526/866554 = 0.9480379
head(overlap)
ewas_EPIC = subset(ewas, CpG %in% overlap) # 42
ewas_EPIC_not_na = subset(ewas_EPIC, !is.na(TWIST)) #27 found by TWIST, 42 - 27 failed QC
ewas_EPIC_not_na = subset(ewas_EPIC_not_na, TWIST < t) # 3
# #              CpG Chromosome  Position         EPIC        TWIST          ONT
# 135674  cg20295214          1 206114536 6.867253e-07 9.676653e-01           NA
# 305247  cg13184736          1  67833726 4.648216e-06 7.580417e-04           NA
# 305255  cg25189904          1  67833810 4.168203e-09 5.556967e-02 4.007738e-05
# 326632  cg09935388          1  92482031 1.568503e-07           NA           NA
# 674072  cg21611682         11  68370801 2.311804e-06 1.001309e-02           NA
# 710749  cg05533761         11  86726911 1.592111e-06           NA           NA
# 710897  cg14391737         11  86802387 7.083859e-16 8.630101e-06           NA
# 710903  cg10711136         11  86804011 1.845844e-09           NA           NA
# 710904  cg00475490         11  86806068 2.245651e-08 1.451271e-02           NA
# 847757  cg09010260         12  54448700 3.433666e-06           NA 1.540150e-03
# 1077441 cg25001882         14  78152734 9.865219e-07           NA 3.023418e-03
# 1122733 cg12586535         15  38822598 4.397695e-06 9.836381e-01           NA

ewas_ONT_not_na_EPIC = subset(suggestive_EPIC, !is.na(ONT)) #21 passed
ewas_ONT_not_na_EPIC = subset(ewas_ONT_not_na_EPIC, ONT < t) #1
dim(ewas_ONT_not_na_EPIC)
#                CpG Chromosome  Position         EPIC        TWIST          ONT
# 305255  cg25189904          1  67833810 4.168203e-09 5.556967e-02 4.007738e-05
# 847757  cg09010260         12  54448700 3.433666e-06           NA 1.540150e-03
# 1077441 cg25001882         14  78152734 9.865219e-07           NA 3.023418e-03
# 1175998 cg18110140         15  75058039 4.215529e-10           NA 1.894146e-03
# 2017834 cg01940273          2 232420224 7.713526e-15 1.014099e-05 6.304819e-03
# 2868082 cg05575921          5    373263 7.958089e-24 2.465543e-11 5.676877e-13
# 2870344 cg04551776          5    393251 2.796015e-06 7.541795e-02 4.326306e-02
# 2870557 cg25648203          5    395329 7.947573e-09 4.555634e-05 1.110545e-03
# 2936211 cg14466441          6  11391960 7.518319e-09           NA 2.055489e-02

# szukamy tego co jest na TWIScie, czego nie zna ONT
overlap = intersect(suggestive_TWIST$CpG, probes)
overlap_TWIST_EPIC = subset(suggestive_TWIST, CpG %in% overlap) # 10
ewas_EPIC_not_na = subset(suggestive_TWIST, !is.na(EPIC)) # 10
ewas_EPIC_not_na = subset(ewas_EPIC_not_na, EPIC < t) # 3

ewas_TWIST_ONT = subset(suggestive_TWIST, !is.na(ONT)) # 71
ewas_TWIST_ONT = subset(ewas_TWIST_ONT, ONT < t) # 12


# szukamy tego co jest na ONT, czego nie zna EPIC

ewas_ONT_not_na_EPIC = subset(suggestive_ONT, !is.na(EPIC)) # 3
ewas_ONT_not_na_EPIC = subset(ewas_ONT_not_na_EPIC, EPIC < t) # 1

ewas_ONT_not_na_TWIST = subset(suggestive_ONT, !is.na(TWIST)) # 17
ewas_ONT_not_na_TWIST = subset(ewas_ONT_not_na_TWIST, TWIST < t) # 12

# total overlap

suggestive_EPIC_TWIST = subset(ewas, TWIST < t & EPIC < t) # 3
suggestive_ONT_EPIC = subset(ewas, ONT < t & EPIC < t) # 1
suggestive_ONT_TWIST = subset(ewas, TWIST < t & ONT < t) # 12

suggestive_EPIC_TWIST_ONT = subset(ewas, TWIST < t & ONT < t & EPIC < t) # 1
