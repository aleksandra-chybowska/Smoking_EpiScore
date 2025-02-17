######################################
## Selection of TWIST samples ########
######################################

## set working directory 
setwd("<local_home_dir>/Desktop/updated_smoking_jan_2019")

## Load requisite libraries 
library(tidyverse)
library(MatchIt)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("lumi")
library(lumi)

### Goal is to select 24 smokers (cases) and 24 never-smokers (controls)
### Smokers to be matched on similar pack years profile 
### Controls to match cases in terms of age and sex 

## read in smoking status information 
smk=read.csv("ever_smoke.csv") 


## Recode smoking status accordingly 
smk[which(smk$ever_smoke==1),"ever_smoke"]<-"Current"
smk[which(smk$ever_smoke==2),"ever_smoke"]<-"Former"
smk[which(smk$ever_smoke==3),"ever_smoke"]<-"Former"
smk[which(smk$ever_smoke==4),"ever_smoke"]<-"Never"
smk[which(smk$ever_smoke==5),"ever_smoke"]<- NA


## read in pack years information 
pck=read.csv("pack_years.csv") 


## Subset datasets to just those with w3 (unrelated DNAm data - will use to form TWIST sub-cohort)
# read in data
w3=read.csv("samplesheet.final.csv")
# subset step 
smk=smk[which(smk$Sample_Name %in% w3$Sample_Name),]
pck=pck[which(pck$Sample_Name %in% w3$Sample_Name),]

#nrow(smk) #4450
#nrow(pck) #4450

## Plot distribution of pack years in current smokers (guide selection of cases)
tot=merge(smk,pck,by="Sample_Name")
curr=tot[tot$ever_smoke %in% "Current",]
plot(hist(curr$pack_years))

## log(pack_years + 1) to get normal distribution
curr$pack_years_transform <- log(curr$pack_years+1)
plot(hist(curr$pack_years_transform))

# # Create the histogram with overlaid data points
# ggplot(curr, aes(x = pack_years_transform)) +
#   geom_histogram(binwidth = 0.5, fill = "#155da2", color = "black", alpha = 0.7) +  # Histogram
#   geom_jitter(aes(y = 0), width = 0, height = 0.2, color = "red", size = 3) +  # Overlay data points
#   labs(
#     x = "Transformed Pack Years",
#     y = "Count",
#     title = "Histogram of Transformed Pack Years with Overlaid Data Points"
#   ) +
#   theme_minimal()
# 
# 

## Select those with large pack years consumption 
length(which(curr$pack_years_transform > 4 & curr$pack_years_transform < 4.5)) ##40 
sub=curr[which(curr$pack_years_transform > 4 & curr$pack_years_transform < 4.5),]
range(sub$pack_years) 

## Exclude outlier - learned after AHRR sanity check but got included in original 24 cases 
sub=sub[-which(sub$Sample_Name %in% "155809"),]


## Randomly sample 24 individuals from group of smokers + high smoking consumption 
## we want approximate sex balance as well as representative range of ages 

## check sex and age distributions 
sel=merge(sub,w3, by= "Sample_Name")

table(sel$sex)
range(sel$age)

## Sample 24 from case group
set.seed(12345)
cases <- sel %>% 
  group_by(sex) %>% 
  sample_n(12, replace = F) 
table(cases$sex)

cases<-as.data.frame(cases)
range(cases$age)

## Select age and sex-matched controls 
cases$Case = 1

## get controls 
nev=tot[tot$ever_smoke %in% "Never",]
nev=merge(nev,w3,by="Sample_Name")
nev$pack_years_transform <- log(nev$pack_years+1)
nev$Case = 0

## combine 
group=rbind(cases,nev)

## Matching step
set.seed(1.2345)
m=matchit(Case ~ age + sex, group, method = "nearest", ratio = 1)
m1 = as.data.frame(m$match.matrix)

m1$Cases <- row.names(m1)
names(m1)[1] <- "Controls"

## pull out cases and matched controls 
ids=m1$Cases
cs = group[match(ids, row.names(group)),]
ids1=m1$Controls
cont = group[match(ids1, row.names(group)),]

nrow(cs)
nrow(cont)

## Assign pairs
cs$Pair=paste0("case",1:24)
cont$Pair=paste0("cont",1:24)

cs$grp=1:24
cont$grp=1:24


## check characteristics 
mean(cont$age,na.rm = T)
mean(cs$age,na.rm = T)
boxplot(cont$age, cs$age)
boxplot(cont$pack_years, cs$pack_years)

## Final group of 48 
test_control = subset(group, group$Case == 0)
test_control$age_diff_55 = 55 - test_control$age



age_diffs <- abs(cs$age - cont$age)

# Combine into a data frame to display the pair differences
pair_diffs <- data.frame(
  Control = cont$Pair,
  Control_Id = cont$Sample_Name,
  Control_Age = cont$age,
  Control_Sex = cont$sex,
  Case = cs$Pair,
  Case_Id = cs$Sample_Name,
  Case_Age = cs$age,
  Case_Sex = cs$sex,
  Age_Difference = age_diffs
)
pair_diffs_when_method_nearest = pair_diffs
write.csv(pair_diffs_when_method_nearest, "~/Desktop/pair_diffs_when_method_nearest.csv", row.names = F)

# full method does not find any matches

pair_diffs_when_method_optimal = pair_diffs
write.csv(pair_diffs_when_method_optimal, "~/Desktop/pair_diffs_when_method_optimal.csv", row.names = F)


tally=rbind(cs,cont)
tally=tally[order(tally$grp),]
# old=read.csv("../Rob/Other/Sample_Selections/TWIST/Samples.csv")
write.csv(tally, "../Rob/Other/Sample_Selections/TWIST/Samples_matched_28072022.csv", row.names =F)

set.seed(2345)
samps=sample(unique(tally$grp),replace=F)
tmp1=samps[1:4]
tmp2=samps[5:8]
tmp3=samps[9:12]
tmp4=samps[13:16]
tmp5=samps[17:20]
tmp6=samps[21:24]

## make swap for better balance
tmp3[2] <-21
tmp5[1] <- 9

grp1=tally[which(tally$grp %in% tmp1),]
grp1$plate=1
grp2=tally[which(tally$grp %in% tmp2),]
grp2$plate=2
grp3=tally[which(tally$grp %in% tmp3),]
grp3$plate=3
grp4=tally[which(tally$grp %in% tmp4),]
grp4$plate=4
grp5=tally[which(tally$grp %in% tmp5),]
grp5$plate=5
grp6=tally[which(tally$grp %in% tmp6),]
grp6$plate=6
tally1=rbind(rbind(rbind(rbind(rbind(grp1,grp2),grp3),grp4),grp5),grp6)

table(tally1$plate, tally1$sex)
with(tally1, tapply(age, plate, range))

## Make new pair and groups for cleanliness 
tally1$pair_new <- paste0(rep(c("case","cont")),rep(1:24,each=2))
tally1$grp_new <- rep(1:24,each=2)
  
## Read in Archie's file 
s1=read.csv("../Rob/Other/Sample_Selections/TWIST/archie_samples.csv")
s2=merge(s1,tally1[,c("Sample_Name", "pair_new","grp_new","plate")],by.x="id",by.y="Sample_Name")
s2=s2[order(s2$grp_new),]
s2$grp_new=NULL
names(s2)[7] <-"pair"
write.csv(s2, "../Rob/Other/Sample_Selections/TWIST/paired_samples_plate_design_28072022.csv",row.names=F)

## Plots of covariates (are they balanced?) 
setwd("../Rob/TWIST/")


## extract variables to be plotted 
var=c("age","Bcell","CD4T", "CD8T","Gran","NK","Mono")

## set up list to store plots 
plot_list = list()

## loop through variables to make boxplots 
for(i in 1:length(var)){
p=ggplot(data=tally,aes(x=factor(Case),y=tally[,var[[i]]],colour=factor(Case)))+geom_boxplot() + xlab("Group") + ylab(var[[i]]) + scale_x_discrete(labels = c("Controls", "Cases")) + scale_colour_discrete(name = "Group", labels = c("Controls", "Cases")) + theme(legend.title.align = 0.5) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
print(p)
plot_list[[i]] <- p 
} 

## write out all plots in one pdf 
pdf("case_control_plots.pdf")
for(i in 1:length(var)){ 
print(plot_list[[i]])  
} 
dev.off()



###################################
#### AHRR PROBE SANITY CHECK ######
###################################
## Read in AHRR data 
ahrr = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Other/Completed_Collabs/Protein Proxies/AHRR_Probe.csv")
names(ahrr)[2] <- "Sample_Sentrix_ID"

## convert m vals to beta vals 
ahrr$AHRR_beta = m2beta(ahrr$AHRR)


## merge with 48 samples 
ahrr=merge(tally,ahrr,by="Sample_Sentrix_ID")

## Plot AHRR probe by case/control status 
pdf("AHRR_check_outlier_removed.pdf")
p=ggplot(data=ahrr,aes(x=factor(Case),y=AHRR_beta,colour=factor(Case)))+geom_boxplot() + xlab("Group") + ylab("% methylation at cg05575921-AHRR") + scale_x_discrete(labels = c("Controls", "Cases")) + scale_colour_discrete(name = "Group", labels = c("Controls", "Cases")) + theme(legend.title.align = 0.5) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
print(p)
dev.off()


## Find outlier in cases: was 155809    
#ahrr.case=ahrr[ahrr$Case %in% 1,]
#ahrr.case[which.max(ahrr.case$AHRR_beta),]
