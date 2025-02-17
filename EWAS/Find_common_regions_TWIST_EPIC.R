################### Read in annotations ######################

if (examine_regions == "TWIST") {

	bed = readRDS(paste0(home, "/TWIST/bed_new_id.RDS"))

	cpgs = as.data.frame(cpgs)
	split = strsplit(cpgs$cpgs, "-", fixed = TRUE)
	cpgs$chr = sapply(split, "[", 1)
	cpgs$start = as.numeric(sapply(split, "[", 2))
	cpgs$stop = as.numeric(sapply(split, "[", 3))

	#                       cpgs   chr    start     stop
	# 1   chr1-21933857-21933858  chr1 21933857 21933858
	# 2   chr1-92481246-92481247  chr1 92481246 92481247
	# 3   chr1-92481318-92481319  chr1 92481318 92481319
	# 4  chr11-68379417-68379418 chr11 68379417 68379418
	# 5  chr11-68379447-68379448 chr11 68379447 68379448
	# 6  chr11-86632814-86632815 chr11 86632814 86632815
	# 7  chr11-86802387-86802388 chr11 86802387 86802388
	# 8  chr13-78594076-78594077 chr13 78594076 78594077
	# 9  chr14-50944025-50944026 chr14 50944025 50944026
	# 10 chr14-89027429-89027430 chr14 89027429 89027430

	colnames(bed) = c("chr", "start", "stop", "desc", "new_id")

	regions = data.frame(
		chr = character(),
		start = numeric(), 
		stop = numeric(),
		start_cpg = numeric(), 
		stop_cpg = numeric(),
		desc = character(), 
		new_id = character() 
	)

	for(i in 1:nrow(cpgs)) {
		site = cpgs[i, ]
		print(site)
		interesting = subset(bed, start <= site$start & stop >= site$stop & chr == site$chr)
		if (nrow(interesting) == 0) {
			interesting[1, ] = c(0, 0, 0, 0, 0, 0, 0)
		} 
		interesting$start_cpg = site$start
		interesting$stop_cpg = site$stop
		interesting$chr = site$chr
		regions = rbind(regions, interesting)
	}

	write.csv(regions, paste0(home, "TWIST/Manual_EWAS/interesting_regions.csv"))

}

	################ second attempt ###############

if (examine_regions == "EPIC") {

	anno = readRDS("<filespace_marioni_group_dir>/Daniel/EPIC_AnnotationObject_df.rds")
	anno = anno[c(1,2,3,18,22,30,31,32,38,40]) # only interesting cols
	# cpg coords don't match, maybe wrong version of the reference genome

	EPICvsTWIST = read.csv(paste0(home, "/TWIST/Manual_EWAS/EPIC.hg38.manifest.tsv"), sep='\t')
	# source: https://zwdzwd.github.io/InfiniumAnnotation (Basic manifest with mapping information, EPIC )
	regions_EPIC = data.frame(
		CpG_chrm = character(),
		CpG_beg = numeric(), 
		CpG_end = numeric(),
		start_cpg = numeric(), 
		stop_cpg = numeric(),
		Probe_ID = character(), 
		new_id = character() 
	)

	for(i in 1:nrow(cpgs)) {
		site = cpgs[i, ]
		print(site)
		interesting = subset(EPICvsTWIST, CpG_beg <= site$start & CpG_end >= site$stop & CpG_chrm == site$chr)
		if (nrow(interesting) == 0) {
			interesting[1, ] = c(0, 0, 0, 0, 0, 0, 0)
		} 
		interesting$start_cpg = site$start
		interesting$stop_cpg = site$stop
		interesting$CpG_chrm = site$chr
		regions_EPIC = rbind(regions_EPIC, interesting)
	}

	write.csv(regions_EPIC, paste0(home, "TWIST/Manual_EWAS/interesting_regions_EPIC.csv"))
}