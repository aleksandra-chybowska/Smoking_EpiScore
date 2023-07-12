#read in bed
link = "https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed"
bed = as.data.frame(read.table(link, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed$new_id = paste0("twist", rownames(bed))
saveRDS(bed, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/bed_new_id.RDS')

out <- strsplit(as.character(bed$V4), ',') 
do.call(rbind, out)
require(data.table)

dt = rbindlist(
  lapply(out, function(x) data.table(t(x))),
  fill = TRUE
)

blended = cbind(bed, dt)
saveRDS(blended, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/bed_exploded.RDS')
