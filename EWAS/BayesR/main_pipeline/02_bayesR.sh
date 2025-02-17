#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Ola, Elena, Rob

############################################################################
# Script to run BayesRR
############################################################################

# use screen -S W3
cd <cluster_home_dir>/Smoking/BayesRR/results/runs/wave3_age_sex/methylation/

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/multiomics/pack_years_W3_4375_agesex_meth_only_processed.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/pack_years_wave3_genetic_age_sex_4375_not_residualized.csvphen \
--analysis-type preprocess \
--thread 24 \
--thread-spawned 24 \
--marker-cache \
--seed 1 > 20230428_pack_years_w3_preprocess.txt


<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/multiomics/pack_years_W3_4375_agesex_meth_only_processed.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/pack_years_wave3_genetic_age_sex_4375_not_residualized.csvphen \
--analysis-type ppbayes \
--chain-length 20000 \
--burn-in 10000 \
--thin 10 \
--S "0.01,0.1,1.0" \
--mcmc-samples pack_years_W3_4375_agesex_meth_only.csv \
--thread 12 \
--thread-spawned 12 \
--marker-cache \
--seed 1 > 20230428_pack_years_w1_w3_9327.txt

## with covariates - wbc

# use screen -S W3
cd <cluster_home_dir>/Smoking/BayesRR/results/runs/white_blood_cells_complete/

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--fixed_effects <cluster_home_dir>/Smoking/BayesRR/data/basic_pack_years_cov.csv \
--fixedEffectNumber 5 \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/methylation_complete/GS20k_allchrom.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/pack_years_complete_17865_residualized.csvphen \
--analysis-type preprocess \
--thread 24 \
--thread-spawned 24 \
--marker-cache \
--seed 3 > 20240819_pack_years_complete_preprocess.txt

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/methylation_complete/GS20k_allchrom.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/pack_years_complete_17865_residualized.csvphen \
--fixed_effects <cluster_home_dir>/Smoking/BayesRR/data/basic_pack_years_cov.csv \
--fixedEffectNumber 5 \
--analysis-type ppbayes \
--chain-length 10000 \
--burn-in 5000 \
--thin 5 \
--S "0.001,0.01,0.1,1.0" \
--mcmc-samples pack_years_17865_complete.csv \
--thread 12 \
--thread-spawned 12 \
--marker-cache \
--seed 4 > 20240819_pack_years_complete.txt
