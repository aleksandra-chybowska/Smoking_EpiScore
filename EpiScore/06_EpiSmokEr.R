library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)
library(EpiSmokEr)

epismoker_cols = c( "cg00006626", "cg03133799", "cg06442199", "cg10006428", "cg14156792", "cg17453416", "cg24629356", "cg26048448", "cg00075467", "cg03936870", "cg05951221", "cg07721625", "cg10362869", "cg17535283", "cg21450627", "cg22587600", "cg02725398", "cg21594961", "cg00379467", "cg08786370", "cg26086649", "cg06677021", "cg08663909", "cg14727987", "cg23576855", "cg27016106", "cg00846554", "cg01555614", "cg03464017", "cg06126421", "cg10525394", "cg15064086", "cg18315060", "cg19643109", "cg19627238", "cg26169299", "cg02243946", "cg03102898", "cg26558023", "cg27650870", "cg13910813", "cg06715410", "cg12589188", "cg03245590", "cg05323345", "cg06120313", "cg09068031", "cg05979241", "cg10609256", "cg10957001", "cg16117605", "cg18161956", "cg19091257", "cg21733098", "cg07499182", "cg23126342", "cg23942311", "cg26103168", "cg18106898", "cg01273991", "cg22947000", "cg00066239", "cg02806012", "cg12438330", "cg13619177", "cg18268547", "cg19572487", "cg27252019", "cg27594073", "cg15129815", "cg06597652", "cg08955995", "cg09298273", "cg20618441", "cg22316634", "cg22331349", "cg13791092", "cg26029902", "cg01103827", "cg01080924", "cg10531774", "cg13451356", "cg21566642", "cg24079702", "cg25242471", "cg13944838", "cg16867657", "cg17619755", "cg25437304", "cg26694437", "cg07179816", "cg13796381", "cg02431260", "cg13771313", "cg04545296", "cg16775095", "cg23766254", "cg24803719", "cg16113156", "cg23513183", "cg18369516", "cg25221984", "cg06644428", "cg09173768", "cg03847932", "cg05575921", "cg20738735", "cg24854181", "cg20839206", "cg05425699", "cg00639837", "cg05293490", "cg18877361", "cg13626582", "cg02267536", "cg06279276", "cg06394109", "cg27436995", "cg00627029", "cg00593900", "cg16702083")

results <- "<cluster_home_dir>/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final/"
pheno <- read.csv(paste0(results, "lbc36_predictions_pack_years_new_target.tsv"), sep = '\t')
pheno$gender <- pheno$sex
pheno$sex <- ifelse(pheno$sex == "M", 1, 2)
rownames(pheno) = pheno$basename

lbc_target <- readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/targets_3489_bloodonly.rds")
rownames(lbc_target) <- lbc_target$Basename

lbc_mvals <- readRDS("<filespace_marioni_group_dir>/LBC/LBC_methylation/LBC_betas_3489_bloodonly.rds")
# for testing
# weights = read.csv(paste0(results,"weights_W1_W3_W4_1256_methylation_adjusted_set_correct_target.csv"))
# weights = weights[-1, ]
# cols = c(epismoker_cols, rownames(lbc_mvals)) # 1376
# cols = unique(cols) # 1371
# lbc_mvals <- subset(lbc_mvals, rownames(lbc_mvals) %in% cols)
# lbc_mvals <- t(lbc_mvals)
# saveRDS(lbc_mvals, '<cluster_home_dir>/Smoking/Elnet_EpiScore/data/lbc_mvals_EpiScore_testing.RDS')
lbc_mvals <- t(lbc_mvals)

lbc_target_36 <- lbc_target[lbc_target["cohort"] == "LBC36",] # 2797, waves 1, 2, 3, and 4
lbc_mvals_36 <- lbc_mvals[rownames(lbc_target_36),]

beta_data = t(lbc_mvals_36)

cols <- intersect(pheno$basename, colnames(beta_data))
beta_data <- beta_data[,cols]
#> identical(colnames(beta_data), pheno$basename)
# [1] TRUE
rownames(pheno) <- pheno$basename

#Dataset has 120 of 121 CpGs required for smoking status estimation.
result <- epismoker(dataset=beta_data, samplesheet = pheno, method = "SSc")
#Dataset has 120 of 121 CpGs required for smoking status estimation.

write.csv(result, "<cluster_home_dir>/Smoking/BayesR_EpiScore/data/sensitivity/EpiSmokEr.csv")