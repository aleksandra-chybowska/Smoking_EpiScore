library(dplyr)
library(readxl)
library(janitor)
library(descr)
library(ggplot2)
library(ggpubr)
box::use(./assign_vars[...])

transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

wd = "/Volumes/marioni-lab/Ola/Lab/Smoking/plots/"

# Glucose, cholesterol, sodium, potassium, urea, creatinine
biochemistry = read.csv('/Volumes/marioni-lab/Generation_Scotland_data/clinical/biochemistry.csv');

smoking_PCQ = read.csv('/Volumes/marioni-lab/Generation_Scotland_data_Sep2021/PCQ/smoking.csv');
smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==1 & is.na(smoking_PCQ$cigs_day), smoking_PCQ$packs_day * 20, smoking_PCQ$cigs_day)
smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==2, smoking_PCQ$packs_day * 20, smoking_PCQ$cigarettes_per_day) # quit within last year
smoking_PCQ$cigarettes_per_day = ifelse(smoking_PCQ$ever_smoke==3, round_half_up((smoking_PCQ$packs_day * 20) / 2.5), smoking_PCQ$cigarettes_per_day) #23459 # quit more than 1 ago

smoking_PCQ = subset(smoking_PCQ, select=-c(age, started, years))
smoking_PCQ = subset(smoking_PCQ, !is.na(ever_smoke)) # 23192
smoking_PCQ = subset(smoking_PCQ, !is.na(cigs_day)) # 22930
smoking_PCQ = subset(smoking_PCQ, !is.na(pack_years)) #22927

target = readRDS("/Volumes/marioni-lab/Ola/Lab/Test_sets/GS20k_Targets.rds")

df = merge(smoking_PCQ, target, by.x="id", by.y="Sample_Name") # 17594
summary(df)
# 6772 NAs in stopped, other cols don't contain NAs
par(mfrow=c(2,3))
hist(df$pack_years, col="steelblue")
hist(df$cigarettes_per_day, col="steelblue")
hist(df$cigs_day, col="steelblue")

hist(outlierTrim(log(df$pack_years+1)), col="steelblue")
hist(outlierTrim(log(df$cigarettes_per_day+1)), col="steelblue")
hist(outlierTrim(log(df$cigs_day+1)), col="steelblue")

df$pack_years = outlierTrim(log(df$pack_years+1))
df$cigarettes_per_day = outlierTrim(log(df$cigarettes_per_day+1))
df$cigs_day = outlierTrim(log(df$cigs_day+1))
summary(df)

smokers_plot = as.data.frame(freq(df$status, y.axis = "percent",  col="steelblue", ylab="%"))
smokers_plot = smokers_plot[-4,]

# total = ggplot(smokers_plot, aes(x=rownames(smokers_plot), y=Percent)) + 
#   geom_bar(stat = "identity", fill="steelblue") +
#   xlab("")
# total
# ggsave(plot=total, filename=paste0(wd, "total.pdf"))

ex_smokers = subset(df, status == "ex-smoker") #5265
quit_more_than_5ya = subset(ex_smokers, stopped > 5) #2284
stopped_by_interval = as.data.frame(freq(ex_smokers$stopped, y.axis = "percent", col="steelblue", ylab="%"))
stopped_by_interval = stopped_by_interval[-b,] ######

ex_smokers$quit_more_than_5 = ifelse(ex_smokers$stopped > 5, 6, ex_smokers$stopped)


plot = ggplot(df, aes(x=cigarettes_per_day, y=pack_years)) +
  geom_point() +
  geom_density_2d() +
  geom_smooth(method=lm, colour='red') +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = 0) + 
  xlab("Cigarettes per day corrected") +
  ylab("Pack years")
ggsave(plot=plot, filename=paste0(wd, "cigs_per_day_corrected_vs_pack_years.pdf"))

plot = ggplot(df, aes(x=cigs_day, y=pack_years)) +
  geom_point() +
  geom_density_2d() +
  geom_smooth(method=lm, colour='red') +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = 0) + 
  xlab("Cigarettes per day") +
  ylab("Pack years")
ggsave(plot=plot, filename=paste0(wd, "cigs_per_day_vs_pack_years.pdf"))

plot = ggplot(df, aes(x=cigarettes_per_day, y=cigs_day)) +
  geom_point() +
  geom_density_2d() +
  geom_smooth(method=lm, colour='red') +
  stat_cor(aes(label = ..r.label..), method = "spearman", cor.coef.name = "r", size = 7, label.x = 0) + 
  xlab("Cigarettes per day") +
  ylab("Cigs per day (uncorrected)")
ggsave(plot=plot, filename=paste0(wd, "cigs_per_day_vs_corrected_cigs_day.pdf"))


