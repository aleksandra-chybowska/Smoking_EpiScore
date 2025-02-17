library(tidyr)
library(dplyr)
library(readxl)

### Prep the longer file

setwd("<datastore_home_dir>/Lab/Smoking/ONT")
# file = read_xlsx("ONT_final_summary.xlsx")
# ds = separate_longer_delim(file, cols = "Native Barcode(s)", delim = "_")
# write.csv(ds, "ONT_final_summary_longer.csv", row.names = F)

### Read in the final file

final_ds = read.csv("Final_files.csv")
annots = read.csv("2022-07-28-Twist_Sample_Info.csv")
annots = annots[c("Sample_Name", "resid")]

complete = merge(final_ds, annots, by.x="Sample.ID", by.y="resid")
complete$GS.Id = complete$Sample_Name
complete = complete[-12]


## test if the DS looks ok

expected_barcode = paste0("barcode", sprintf("%02d", complete$Native.Barcode))
# path = complete$File.Path
# for (i in 1:length(expected_barcode)) {
#   print(grepl(pattern = expected_barcode[i], x = path[i], ignore.case = TRUE))
# }
print(mapply(grepl, pattern = expected_barcode, x = complete$File.Path, 
             MoreArgs = list(ignore.case = TRUE)))
table(duplicated(complete$GS.Id))
write.csv(complete, "2024-10-04-ONT_Sample_Info.csv", row.names = F)


