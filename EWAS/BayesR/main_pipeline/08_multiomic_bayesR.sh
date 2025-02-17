#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Ola, Elena, Rob

############################################################################
# Script to run BayesRR
############################################################################

# use screen -S W3
# screen -S multiomics_20k
cd <cluster_home_dir>/Smoking/BayesRR/results/runs/w1_w3_genetics_age_sex/multiomics/

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/multiomics_W1_W3/pack_years_W3_W4_allchrom_processed.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/pack_years_wave1_wave3_genetic_age_sex_9316_not_residualized.csvphen \
--group <cluster_home_dir>/Smoking/BayesRR/data/multiomics_W1_W3/pack_years_cpg_groups.txt \
--analysis-type preprocess \
--thread 1 \
--S "0.0001,0.001,0.01;0.01,0.1,1.0" \
--thread-spawned 1 \
--marker-cache \
--seed 1 > 20230429_pack_years_W3_W1_preprocess.txt

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/multiomics_W1_W3/pack_years_W3_W4_allchrom_processed.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/pack_years_wave1_wave3_genetic_age_sex_9316_not_residualized.csvphen \
--group <cluster_home_dir>/Smoking/BayesRR/data/multiomics_W1_W3/pack_years_cpg_groups.txt \
--analysis-type ppbayes \
--chain-length 20000 \
--burn-in 10000 \
--thin 10 \
--S "0.0001,0.001,0.01;0.01,0.1,1.0" \
--mcmc-samples pack_years_W1_W3_9316_age_sex_genetic_cpg_seed_test.csv \
--thread 12 \
--thread-spawned 12 \
--marker-cache \
--seed 234234 > 20230505_pack_years_W1_W3_9316_age_sex_genetic_cpg_large_seed.txt
