link = "https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed"
bed <- as.data.frame(read.table(link, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))

out <- strsplit(as.character(bed$V4),',') 
do.call(rbind, out)
require(data.table)

dt = rbindlist(
  lapply(out, function(x) data.table(t(x))),
  fill = TRUE
)

