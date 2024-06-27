#Translating 

#EUR stratified without UK Biobank

./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified_without_UKB/EUR_stratified_without_UKB/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt \
--out /home/shirin/data/EUR_stratified_without_UKB/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified_without_UKB/EUR_stratified_without_UKB/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt \
--out /home/shirin/data/EUR_stratified_without_UKB/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified_without_UKB/EUR_stratified_without_UKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt \
--out /home/shirin/data/EUR_stratified_without_UKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist


#=====
#EUR stratified
./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified/EUR_stratified/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out /home/shirin/data/EUR_stratified/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified/EUR_stratified/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out /home/shirin/data/EUR_stratified/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified/EUR_stratified/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out /home/shirin/data/EUR_stratified/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/EUR_stratified/EUR_stratified/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt \
--out /home/shirin/data/EUR_stratified/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

##weird alleles?
./munge_sumstats.py \
--sumstats /home/shirin/data/grimage_translated.outfile.fastGWA.tab \
--out /home/shirin/data/grimage_PACKYRS \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/grimage_translated_alleles_swapped.outfile.fastGWA.tab \
--out /home/shirin/data/grimage_translated_alleles_swapped \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

./munge_sumstats.py \
--sumstats /home/shirin/data/ErzurumluogluAM_30617275_discovery-stage_meta-analysis_PackYears.txt \
--out /home/shirin/data/largest \
--merge-alleles /home/shirin/data/eur_w_ld_chr/w_hm3.snplist

# Largest

./ldsc.py \
--rg /home/shirin/data/EpiSmokEr.sumstats.gz,/home/shirin/data/largest.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/epismoker_largest

# Grimage vs EUR stratified without UK Biobank
./ldsc.py \
--rg /home/shirin/data/grimage_PACKYRS.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_PACKYRS_EUR_without_UKB_CigDay

./ldsc.py \
--rg /home/shirin/data/grimage_translated_alleles_swapped.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_translated_alleles_swapped_PACKYRS_EUR_without_UKB_CigDay

# negative correlation
./ldsc.py \
--rg /home/shirin/data/grimage_translated_alleles_swapped.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_PACKYRS_EUR_without_UKB_AgeSmk

./ldsc.py \
--rg /home/shirin/data/grimage_PACKYRS.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_PACKYRS_EUR_without_UKB_SmkCes

./ldsc.py \
--rg /home/shirin/data/grimage_PACKYRS.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_PACKYRS_EUR_without_UKB_SmkInit


# bulk - this is great!

./ldsc.py \
--rg /home/shirin/data/grimage_PACKYRS.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_PACKYRS_EUR_without_UKB_all

./ldsc.py \
--rg /home/shirin/data/grimage_PACKYRS.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/grimage_PACKYRS_EUR_stratified_all

./ldsc.py \
--rg /home/shirin/data/EpiSmokEr.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz,/home/shirin/data/EUR_stratified_without_UKB/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/EpiSmokEr_translated_EUR_without_UKB_all

./ldsc.py \
--rg /home/shirin/data/EpiSmokEr.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz \
--ref-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--w-ld-chr /home/shirin/data/eur_w_ld_chr/ \
--out /home/shirin/data/EpiSmokEr_translated_EUR_stratified_all

/home/shirin/data/EUR_stratified/munged/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz,/home/shirin/data/EUR_stratified/munged/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.sumstats.gz