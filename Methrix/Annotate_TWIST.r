annots = readRDS('/Cluster_Filespace/Marioni_Group/Ola/TWIST/Annotations/40_10_CpGs_annotated.RDS')
ewas = readRDS('/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/twist/40_10_EWAS_stats_twist.RDS')
head(ewas)

ewas_significant = subset(ewas, ewas$pval < 1e-5) # 97
annots = annots[c(1:3, 54:63)]
annots = annots %>% unite("cpgs", c("seqnames", "start", "end"), sep="-")
saveRDS(annots, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/Annotations/40_10_CpGs_annotated_min.RDS', compress = F)
df = merge(ewas_significant, annots, by.x = "cpg", by.y="cpgs")
write.csv(df, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/Annotations/40_10_annotated_hits.csv', row.names=F)

############
# Now EPIC

anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
anno <- anno[c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
ewas <- readRDS('/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/epic/40_10_EWAS_stats_epic.RDS')
head(ewas)
ewas_significant = subset(ewas, ewas$pval < 1e-5) #44
epic = merge(ewas_significant, anno, by.x = "cpg", by.y = "Name")
write.csv(epic, '/Cluster_Filespace/Marioni_Group/Ola/TWIST/Manual_EWAS/epic/40_10_annotated_hits.csv', row.names=F)
