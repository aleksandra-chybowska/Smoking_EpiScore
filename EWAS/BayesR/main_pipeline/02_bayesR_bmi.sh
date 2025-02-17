#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Ola, Elena, Rob

export LD_LIBRARY_PATH=/opt/gcc/lib64

############################################################################
# Script to run BayesRR
############################################################################

screen -S All-1
cd <cluster_home_dir>/Smoking/BayesRR/results/runs/bmi

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/methylation/GS20k_allchrom.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/bmi_17833_residualized.csvphen \
--analysis-type preprocess \
--thread 24 \
--thread-spawned 24 \
--marker-cache \
--seed 1 > 20230317_bmi.txt

<filespace_marioni_group_dir>/BayesRRcmd/src/brr \
--data-file <cluster_home_dir>/Smoking/BayesRR/data/methylation/GS20k_allchrom.csv \
--pheno <cluster_home_dir>/Smoking/BayesRR/data/bmi_17833_residualized.csvphen \
--analysis-type ppbayes \
--chain-length 10000 \
--S "0.0001,0.001,0.01" \
--mcmc-samples bmi_17833.csv \
--thread 12 \
--thread-spawned 12 \
--marker-cache \
--seed 1 > 20230317_bmi_ewas.txt

