anno2 = fread("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Cohort_Summary_Stats/EasyQC/Allele_Freq_and_Mapping_Info/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.", stringsAsFactors=F, header=T)
anno = fread("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Cohort_Summary_Stats/EasyQC/Allele_Freq_and_Mapping_Info/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map.gz", stringsAsFactors=F, header=T)

anno$cptid = paste0(anno$chr, ":", anno$pos)
anno_full = merge(anno, anno2, by.x = "cptid", by.y = "cptid")
anno_full$match = paste0(anno_full$chr, "_", anno_full$pos, "_", anno_full$ref, "_", anno_full$alt)
fwrite(anno_full, "/Cluster_Filespace/Marioni_Group/Ola/Smoking/GWAS/anno_full.csv", sep=",", row.names=F, col.names=F, quote=F)