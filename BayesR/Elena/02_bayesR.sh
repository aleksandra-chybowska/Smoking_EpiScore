#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Ola, Elena, Rob

############################################################################
# Script to run BayesRR
############################################################################

# use screen -S W3
cd /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/results/runs/wave3_age_sex/methylation/

/Cluster_Filespace/Marioni_Group/BayesRRcmd/src/brr \
--data-file /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/multiomics/pack_years_W3_4375_agesex_meth_only_processed.csv \
--pheno /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_wave3_genetic_age_sex_4375_not_residualized.csvphen \
--analysis-type preprocess \
--thread 24 \
--thread-spawned 24 \
--marker-cache \
--seed 1 > 20230428_pack_years_w3_preprocess.txt


/Cluster_Filespace/Marioni_Group/BayesRRcmd/src/brr \
--data-file /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/multiomics/pack_years_W3_4375_agesex_meth_only_processed.csv \
--pheno /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/data/pack_years_wave3_genetic_age_sex_4375_not_residualized.csvphen \
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

