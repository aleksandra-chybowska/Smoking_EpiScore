library("ggplot2")
library("stringr")
library("data.table")
library("foreign")

results <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Elnet_EpiScore/results/j_1e-4_pack_years_20k_final/"
path <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/data/"
neuro_path <- "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Neuroimaging/"

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

## Load brain data
####################################################################

pheno <- read.csv(paste0(results, "lbc36_predictions_pack_years.tsv"), sep = '\t')
pheno_w1 <- pheno[pheno$WAVE == 1,] # 882

ph36 <- read.spss(paste0(path,
    "LBC/AleksandraChybowska_BloodAndBrainBasedEWASofSmoking_AC_05APR2023.sav"),
    to.data.frame = TRUE,
    use.value.labels = FALSE)

ph36_min <- ph36[c("lbc36no", "ICV_mm3_wX", "brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2")]
ph36_min$lbc36no <- trimws(ph36_min$lbc36no)

rownames(ph36_min) <- ph36_min$lbc36no
ph36_min <- ph36_min[-c(1)]
ph36_min_ICV <- ph36_min / ph36_min$ICV_mm3_wX
ph36_min$lbc36no <- rownames(ph36_min)
ph36_min_ICV$lbc36no <- rownames(ph36_min_ICV)

neuroimaging <- merge(pheno_w1, ph36_min, by.x="ID", by.y="lbc36no")
norm_neuroimaging <- merge(pheno_w1, ph36_min_ICV, by.x="ID", by.y="lbc36no")

summary(lm(scale(brain_mm3_w2) ~ scale(age) + as.factor(sex) + scale(py_pred) + scale(ICV_mm3_wX), data = neuroimaging))
summary(lm(scale(brain_mm3_w2) ~ scale(age) + as.factor(sex) + scale(pack_years_clean) + scale(ICV_mm3_wX), data = neuroimaging))

saveRDS(neuroimaging, paste0(neuro_path, "neuroimaging_pheno.RDS"), compress = F)
saveRDS(norm_neuroimaging, paste0(neuro_path, "neuroimaging_pheno_ICV_norm.RDS"), compress = F)

## Variable distribution check, remove NAs, scale
####################################################################

neuroimaging <- neuroimaging[c("ID",  "age", "sex", "py_pred", "pack_years_clean",
                      "ICV_mm3_wX", "brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2")]
neuroimaging <- na.omit(neuroimaging) #525

columns = c("brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2", "ICV_mm3_wX", "py_pred", "pack_years_clean")
pdf(paste0(neuro_path, "outcome_distribution_neuroimaging.pdf"))
par(mfrow = c(3, 2))
for(col in columns) {
  hist(neuroimaging[,col], main = col, breaks=20, xlab="")
}
dev.off()

print("Scaling variables, not tranforming them - watch out for wmh.")
neuroimaging$brain_mm3_w2_sc_t = scale(neuroimaging$brain_mm3_w2)
neuroimaging$wmh_mm3_w2_sc_t = scale(transform(neuroimaging$wmh_mm3_w2))
neuroimaging$gm_mm3_w2_sc_t = scale(neuroimaging$gm_mm3_w2)
neuroimaging$nawm_mm3_w2_sc_t = scale(neuroimaging$nawm_mm3_w2)
neuroimaging$ICV_mm3_wX_sc_t = scale(neuroimaging$ICV_mm3_wX)
neuroimaging$py_pred_sc = scale(neuroimaging$py_pred)
neuroimaging$pack_years_clean_sc = scale(neuroimaging$pack_years_clean)

out_cont <- data.frame(Outcome=character(), 
	Predictor=character(),
	n=double(),
	Beta=double(),
	SE=double(),
	P=double(), 
	LCI=double(),
	UCI=double(),
	stringsAsFactors=FALSE)

my_data = neuroimaging
y <- c("brain_mm3_w2_sc_t", "wmh_mm3_w2_sc_t", "gm_mm3_w2_sc_t", "nawm_mm3_w2_sc_t")
x <- c("py_pred_sc", "pack_years_clean_sc") 

count <- 0

for(i in 1:length(y)) {
	for(j in 1:length(x)) {

	count = count + 1
	
	outcome <- y[i]
	predictor <- x[j]
	
	model <- lm(my_data[,y[i]] ~ my_data[,x[j]] + scale(age) + as.factor(sex) + ICV_mm3_wX_sc_t, data=my_data) # to be edited manually for covariates 
	coefs <- summary(model)$coefficients[2,c(1,2,4)] # need to make sure that the coefficient of interest is first in the model
	int <- confint(model)[2,1:2]
	n = nobs(model)

	out_cont[count, ] <- c(outcome, predictor, n, signif(coefs, 2), signif(int, 2)) 

	}
}

# Daniel Excludes people with signs of dementia and uses wave2 for associations
out_cont <- as.data.frame(out_cont)
out_cont$Beta <- as.numeric(out_cont$Beta)
out_cont$SE <- as.numeric(out_cont$SE)
out_cont$LCI <- as.numeric(out_cont$LCI)
out_cont$UCI <- as.numeric(out_cont$UCI)
saveRDS(out_cont, paste0(neuro_path, "trasformed_lm_res.RDS"), compress = F)

## Prep plot
####################################################################
out_cont <- readRDS(paste0(neuro_path, "trasformed_lm_res.RDS"))
out_cont$Outcome = gsub("brain_mm3_w2_sc_t", "Total brain\nvolume", out_cont$Outcome)
out_cont$Outcome = gsub("wmh_mm3_w2_sc_t", "White matter\nhyperintensity\nvolume", out_cont$Outcome)
out_cont$Outcome = gsub("gm_mm3_w2_sc_t", "Gray matter\nvolume", out_cont$Outcome)
out_cont$Outcome = gsub("nawm_mm3_w2_sc_t", "Normal appearing\nwhite matter\nvolume", out_cont$Outcome)
out_cont$Predictor = gsub("py_pred_sc", "EpiScore", out_cont$Predictor)
out_cont$Predictor = gsub("pack_years_clean_sc", "Pack years", out_cont$Predictor)

My_Theme = theme(
  panel.border = element_rect(colour="black",size=1, fill = NA),
  axis.title.x = element_text(size = 16), # controls HR label size 
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 16, face = "bold"),
  legend.text=element_text(size=16),
	    legend.position=c(0.8,0.1),
	    legend.background = element_rect(size=0.5, linetype="solid", colour="black"),
	    legend.title=element_blank(),
  axis.title=element_text(size=16))

stacked = ggplot(out_cont,aes(y=Beta, x=Outcome, group=Predictor, colour=Predictor)) + 
  geom_point(size = 2, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1)+
  ylab("Effect Size (Beta)")+ 
  xlab ("") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5), axis.text.y = element_text(size = 8), legend.position = "right",
        plot.title = element_text(size = 8))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + My_Theme


pdf(paste0(neuro_path, "Forest_plot_neuroimaging_not_transformed_wmh_transformed.pdf"))
stacked
dev.off()