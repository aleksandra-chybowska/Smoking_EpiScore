library("tidyverse")
setwd("U:/Ola/Lab/Smoking/ONT/")
dict = read_csv("2022-07-28-Twist_Sample_Info.csv")
files = read_csv("Paths_and_IDS.csv")

info =  merge(files, dict, by.x="Sample_ID", by.y="resid")
to_repeat = read_csv("CpGsToRepeat.csv", col_names = F)
intersect(to_repeat$X1, info$Sample_ID)



merged = merge(files, dict[c("Sample_Name", "resid")], by.x="Sample_ID", by.y="resid")
merged$new_name = paste0("./", merged$Sample_Name, ".wf_mods.bedmethyl.gz")
write_csv(merged, "renaming_strategy.csv")
