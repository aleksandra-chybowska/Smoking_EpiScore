##################################################################################
### Script to generate DNAm predictor of Smoking in LBC1936 at wave 1 (age 70) ###
##################################################################################

### Look at phenotypes - choose outcome ###
library("foreign")
library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("caret")

setwd("<cluster_home_dir>/Smoking/BayesR_EpiScore/data/")

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

### calculate pack_years based on cigarettes per day ###

cigs <- read.spss(
    "LBC/LBC1936_EWAS_BldBasedChonicLowGradeInflammation_RH_27FEB2023.sav",
    to.data.frame = TRUE,
    use.value.labels = FALSE)

cigs_min <- cigs[c("lbc36no",
                    "smokagestart_w1",
                    "smokagestop_w1",
                    "smoknumcigs_w1")]

ph36 <- read.spss(
    "LBC/AleksandraChybowska_BloodAndBrainBasedEWASofSmoking_AC_05APR2023.sav",
    to.data.frame = TRUE,
    use.value.labels = FALSE)

ph36$lbc36no <- trimws(ph36$lbc36no)

ph36_min <- ph36[c("lbc36no",
                    "sex",
                    "agedays_w1",
                    "smokcat_w1")]

ph <- merge(ph36_min, cigs_min, by.x = "lbc36no", by.y = "lbc36no")
dim(ph)
# 1091

### Check if there are any people who started smoking but stop date is NA - Show to Riccardo
### 0 - never smokers, 1 - ex-smoker, 2 - current smoker

subset(ph, !is.na(smokagestart_w1) & is.na(smokagestop_w1) & smokcat_w1 != 2)
to_remove <- which(!is.na(ph$smokagestart_w1) &
                    is.na(ph$smokagestop_w1) &
                    ph$smokcat_w1 != 2)
ph <- ph[-to_remove, ]
dim(ph)
# 1081

subset(ph, !is.na(smokagestart_w1) & is.na(smokagestop_w1))
index <- which(!is.na(ph$smokagestart_w1) & is.na(ph$smokagestop_w1))
ph[index, "smokagestop_w1"] <- 70

### Number of cigarettes smoked per day / pack size * number of years a person has smoked
ph$pack_years <- (ph$smokagestop_w1 - ph$smokagestart_w1) * ph$smoknumcigs_w1/20
ph[which(ph$smokcat_w1 == 0), "pack_years"] <- 0
#summary(ph)

### Where do NAs come from?
#subset(ph, is.na(pack_years))
### Former smokers who dont know how much they smoked - safe to remove

ph <- subset(ph, !is.na(pack_years))
ph$pack_years_clean <- outlierTrim(log(ph$pack_years + 1))
ph <- subset(ph, !is.na(pack_years_clean))
#dim(ph)
### 1072

ph$ever_smoke <- as.factor(ph$smokcat_w1)
sample_size <- nrow(ph)

pdf("LBC/pheno_min_hist.pdf")
hist(ph$pack_years, breaks = 20)
dev.off()

pdf("LBC/pheno_min_log_hist.pdf")
hist(ph$pack_years_clean, breaks = 20)
dev.off()

out <- paste0("LBC/pheno_min_", sample_size, ".csv")
write.csv(ph, out, row.names = FALSE)
print("Saved!")
print(out)
### Finished prepping pheno