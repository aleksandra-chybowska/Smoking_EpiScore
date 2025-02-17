#Translating 

#EUR stratified without UK Biobank

./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified_without_UKB/EUR_stratified_without_UKB/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt \
--out <local_home_dir>/data/EUR_stratified_without_UKB/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified_without_UKB/EUR_stratified_without_UKB/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt \
--out <local_home_dir>/data/EUR_stratified_without_UKB/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified_without_UKB/EUR_stratified_without_UKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt \
--out <local_home_dir>/data/EUR_stratified_without_UKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist


#=====
#EUR stratified
./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified/EUR_stratified/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out <local_home_dir>/data/EUR_stratified/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified/EUR_stratified/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out <local_home_dir>/data/EUR_stratified/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified/EUR_stratified/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out <local_home_dir>/data/EUR_stratified/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/EUR_stratified/EUR_stratified/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out <local_home_dir>/data/EUR_stratified/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

##weird alleles?
./munge_sumstats.py \
--sumstats <local_home_dir>/data/grimage_translated.outfile.fastGWA.tab \
--out <local_home_dir>/data/grimage_PACKYRS \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/grimage_translated_alleles_swapped.outfile.fastGWA.tab \
--out <local_home_dir>/data/grimage_translated_alleles_swapped \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats <local_home_dir>/data/ErzurumluogluAM_30617275_discovery-stage_meta-analysis_PackYears.txt \
--out <local_home_dir>/data/largest \
--merge-alleles <local_home_dir>/data/eur_w_ld_chr/w_hm3.snplist

# Largest

./ldsc.py \
--rg <local_home_dir>/data/EpiSmokEr.sumstats.gz,<local_home_dir>/data/largest.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/epismoker_largest

# Grimage vs EUR stratified without UK Biobank
./ldsc.py \
--rg <local_home_dir>/data/grimage_PACKYRS.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_PACKYRS_EUR_without_UKB_CigDay

./ldsc.py \
--rg <local_home_dir>/data/grimage_translated_alleles_swapped.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_translated_alleles_swapped_PACKYRS_EUR_without_UKB_CigDay

# negative correlation
./ldsc.py \
--rg <local_home_dir>/data/grimage_translated_alleles_swapped.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_PACKYRS_EUR_without_UKB_AgeSmk

./ldsc.py \
--rg <local_home_dir>/data/grimage_PACKYRS.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_PACKYRS_EUR_without_UKB_SmkCes

./ldsc.py \
--rg <local_home_dir>/data/grimage_PACKYRS.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_PACKYRS_EUR_without_UKB_SmkInit


# bulk - this is great!

./ldsc.py \
--rg <local_home_dir>/data/grimage_PACKYRS.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_PACKYRS_EUR_without_UKB_all

./ldsc.py \
--rg <local_home_dir>/data/grimage_PACKYRS.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/grimage_PACKYRS_EUR_stratified_all

./ldsc.py \
--rg <local_home_dir>/data/EpiSmokEr.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,<local_home_dir>/data/EUR_stratified_without_UKB/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/EpiSmokEr_translated_EUR_without_UKB_all

./ldsc.py \
--rg <local_home_dir>/data/EpiSmokEr.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz \
--ref-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--w-ld-chr <local_home_dir>/data/eur_w_ld_chr/ \
--out <local_home_dir>/data/EpiSmokEr_translated_EUR_stratified_all

<local_home_dir>/data/EUR_stratified/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,<local_home_dir>/data/EUR_stratified/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz