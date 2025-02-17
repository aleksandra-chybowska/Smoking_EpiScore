#cd <filespace_marioni_group_dir>/Anna/Brain_methylation/

### Move to R ###
setwd("<filespace_marioni_group_dir>/Anna/Brain_methylation/")
dat = readRDS("Normalised_Beta_Values_danet.rds")
descriptions = readRDS("descriptions.rds") # According to this file, this patient does not have HC data: LBC360666

dat <- t(dat)
identical(rownames(descriptions), rownames(dat))

dat <- merge(descriptions, dat, by="row.names")

library("foreign")
ph36 <- read.spss("<filespace_marioni_group_dir>/Riccardo/Cognitive_EWAS_GS_LBC/LBC1936_BloodBasedEWASofCognitiveAbilities_RM_06MAY2021.sav", to.data.frame=T, use.value.labels=F)
ph361 <- ph36[,c(1,12)]
names(ph361)[1] <- "LBC_ID"

d <- merge(ph361, dat, by="LBC_ID")

d1 <- d[-grep("_R", d$Sample_Name),] #  69 807184

hc <- d1[grep("HC", d1$Sample_Name),] # smokcat_w1 missing for one person with hc
pfc <- d1[grep("BA46", d1$Sample_Name),]
ba17 <- d1[grep("BA17", d1$Sample_Name),]
ba24 <- d1[grep("BA24", d1$Sample_Name),]
ba2021 <- d1[grep("BA20", d1$Sample_Name),]

record_list = list(
	hc=list(hc$cg05575921, hc$smokcat_w1), 
	pfc=list(pfc$cg05575921, pfc$smokcat_w1),
	ba17=list(ba17$cg05575921, ba17$smokcat_w1),
	ba24=list(ba24$cg05575921, ba24$smokcat_w1),
	ba2021=list(ba2021$cg05575921, ba2021$smokcat_w1)
)
hc_interesting = hc[c("LBC_ID", "smokcat_w1",
	"cg09558034", "cg26381592", "cg10372485", 
	"cg10892497", "cg08545123", "cg05575921")]

saveRDS(hc_interesting, "<cluster_home_dir>/Smoking/Brain_vs_Blood/hc_interesting.R")
saveRDS(record_list, "<cluster_home_dir>/Smoking/Brain_vs_Blood/AHRR_record_list.R", compress=F)

### AHRR plots ###
pdf(file="<cluster_home_dir>/Smoking/Brain_vs_Blood/AHRR.pdf")
par(mfrow=c(2,3))
	plot(hc$cg05575921 ~ hc$smokcat_w1)
	plot(pfc$cg05575921 ~ pfc$smokcat_w1)
	plot(ba17$cg05575921 ~ ba17$smokcat_w1)
	plot(ba24$cg05575921 ~ ba24$smokcat_w1)
	plot(ba2021$cg05575921 ~ ba2021$smokcat_w1)
dev.off()

hc_p <- NA
hc_cg <- NA
for(i in 22:ncol(hc)){
hc_p[i-21] <- summary(lm(hc[,i] ~ hc$smokcat_w1))$coefficients[2,4]
hc_cg[i-21] <- names(hc)[i]
}


pfc_p <- NA
pfc_cg <- NA
for(i in 22:ncol(pfc)){
pfc_p[i-21] <- summary(lm(pfc[,i] ~ pfc$smokcat_w1))$coefficients[2,4]
pfc_cg[i-21] <- names(pfc)[i]
}

ba17_p <- NA
ba17_cg <- NA
for(i in 22:ncol(ba17)){
ba17_p[i-21] <- summary(lm(ba17[,i] ~ ba17$smokcat_w1))$coefficients[2,4]
ba17_cg[i-21] <- names(ba17)[i]
}

ba24_p <- NA
ba24_cg <- NA
for(i in 22:ncol(ba24)){
ba24_p[i-21] <- summary(lm(ba24[,i] ~ ba24$smokcat_w1))$coefficients[2,4]
ba24_cg[i-21] <- names(ba24)[i]
}

ba2021_p <- NA
ba2021_cg <- NA
for(i in 22:ncol(ba2021)){
ba2021_p[i-21] <- summary(lm(ba2021[,i] ~ ba2021$smokcat_w1))$coefficients[2,4]
ba2021_cg[i-21] <- names(ba2021)[i]
}

results_hc = data.frame(
	"hc_p" = hc_p, 
	"hc_cg" = hc_cg
)

results_pfc = data.frame(
	"pfc_p" = pfc_p,
	"pfc_cg" = pfc_cg
)

results_ba17 = data.frame(
	"ba17_p" = ba17_p,
	"ba17_cg" = ba17_cg
)

results_ba24 = data.frame(
	"ba24_p" = ba24_p,
	"ba24_cg" = ba24_cg
)

results_ba2021 = data.frame(
	"ba2021_p" = ba2021_p,
	"ba2021_cg" = ba2021_cg
)
path <- "<cluster_home_dir>/Smoking/Brain_vs_Blood/"
saveRDS(results_hc, paste0(path, "results_hc.RDS"), compress=F)
saveRDS(results_pfc, paste0(path, "results_pfc.RDS"), compress=F)
saveRDS(results_ba17, paste0(path, "results_ba17.RDS"), compress=F)
saveRDS(results_ba24, paste0(path, "results_ba24.RDS"), compress=F)
saveRDS(results_ba2021, paste0(path, "results_ba2021.RDS"), compress=F)

saveRDS(hc, paste0(path, "hc.RDS"), compress=F)
saveRDS(pfc, paste0(path, "pfc.RDS"), compress=F)
saveRDS(ba17, paste0(path, "ba17.RDS"), compress=F)
saveRDS(ba24, paste0(path, "ba24.RDS"), compress=F)
saveRDS(ba2021, paste0(path, "ba2021.RDS"), compress=F)

sig1 <- length(which(hc_p < 1e-5))
sig11 <- which(hc_p < 1e-5)
tmp <- sig11 + 21
cpgs <- names(hc)[tmp]
# "cg09558034" "cg26381592" "cg10372485" "cg10892497" "cg08545123"

cpgs_from_hc = hc[tmp]
cpgs_from_hc$smoking_cat = hc$smokcat_w1
#saveRDS(cpgs_from_hc, "<cluster_home_dir>/Smoking/Plots/HC.RDS")
pdf(file="<cluster_home_dir>/Smoking/Plots/HC.pdf")
par(mfrow=c(2,3))
for(i in 1:sig1){
plot(hc[,sig11[i] + 21] ~ hc$smokcat_w1, main=cpgs[i])

}
dev.off()
sig1 <- length(which(ba2021_p < 1e-5))
sig11 <- which(ba2021_p < 1e-5)
tmp <- sig11 + 21
cpgs <- names(ba2021)[tmp]

par(mfrow=c(1,sig1))
for(i in 1:sig1){
plot(ba2021[,sig11[i] + 21] ~ ba2021$smokcat_w1, main=cpgs[i])
}
