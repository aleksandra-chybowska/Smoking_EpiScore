library(data.table)
library(tidyverse)

setwd("<datastore_home_dir>/Lab/Smoking/Replication")

## ALSPAC
##   time_code     n mean.age sd.age n.female percent.female n.current n.former
##   <chr>     <int>    <dbl>  <dbl>    <int>          <dbl>     <int>    <int>
## 1 F17        1214     17.7  0.372      662           54.5       110      132
## 2 F24         496     24.4  0.762      264           53.2       155       96
## 3 FOM_FOF    1207     50.2  5.38       668           55.3        90      474
## 4 antenatal   968     28.8  4.44       968          100          98      268

extract_summary = function(df, name, current_flag=2, former_flag=1) {
  n = nrow(df)
  mean.age = mean(df$age)
  sd.age = sd(df$age)
  n.female = nrow(df[df$sex=="F"])
  n.male = n - n.female
  percent.female = (n.female / n) * 100
  percent.male = 100 - percent.female
  n.current = nrow(df[df$smokcat_w1 == current_flag])
  n.former = nrow(df[df$smokcat_w1 == former_flag])
  n.never = n - (n.current + n.former)
  
  ret = data.frame(n, mean.age, sd.age, n.female, n.male, percent.female, 
                   percent.male, n.current, n.former, n.never)
  return(ret)
}

GS=fread("pack_years_17865_complete.csv")
GS_smk=fread("GS_smoking_statuses.csv")

GS=merge(GS, GS_smk, by="Sample_Name")
GS=rename(GS, smokcat_w1=new_group)
GS_pheno=extract_summary(GS, "Waves 1, 3 and 4")
GS_pheno$study="GS"
GS_pheno$n.450k=0
GS_pheno$n.EPIC=GS_pheno$n
GS_pheno$time_code="main_GS"

LBC=fread("LBC_phenos_summary.csv") 
LBC_pheno=extract_summary(LBC, "LBC_main")
LBC_pheno$study="LBC1931"
LBC_pheno$n.450k = LBC_pheno$n
LBC_pheno$n.EPIC = 0
LBC_pheno$time_code="Wave 1"

ALSPAC=fread("ALSPAC.demo.tsv")
ALSPAC$study="ALSPAC"
ALSPAC$n.male = ALSPAC$n - ALSPAC$n.female
ALSPAC$percent.male = 100 - ALSPAC$percent.female
ALSPAC$n.never = ALSPAC$n - (ALSPAC$n.current + ALSPAC$n.former)

summary=rbind(GS_pheno, LBC_pheno, ALSPAC)
fwrite(summary, "summary.csv")
