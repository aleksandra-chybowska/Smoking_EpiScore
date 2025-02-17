library("tidyverse")

setwd("<marionettes_data_dir>/nanopore_v_lmurphy2/extracted_bedgraphs")
path <- "."

# List all files and directories recursively
files_and_dirs <- list.files(path = path, recursive = TRUE,
                            full.names = TRUE, pattern = "*.bedmethyl.gz")

# Print the file hierarchy
print(files_and_dirs)
write.csv(files_and_dirs, "File_list.csv")

### Link with ids

setwd("<datastore_home_dir>/Lab/Smoking/ONT/")
dict = read_csv("2022-07-28-Twist_Sample_Info.csv")
files = read_csv("Paths_and_IDS.csv")
merged = merge(files, dict[c("Sample_Name", "resid")], by.x="Sample_ID", by.y="resid")
merged$new_name = paste0("./", merged$Sample_Name, ".wf_mods.bedmethyl.gz")
write_csv(merged, "renaming_strategy.csv")

### Rename files based on renaming strategy
setwd("<marionettes_data_dir>/nanopore_v_lmurphy2/extracted_bedgraphs/")
renaming_strategy = read_csv("renaming_strategy.csv")

# Copy from old files to new files
file.rename(from = renaming_strategy$File, to = renaming_strategy$new_name)
