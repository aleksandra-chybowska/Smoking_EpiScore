cd /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/

../../../BayesRRcmd/src/brr \
--data-file Inputs/Methylation/pack_years_17833.csv \
--pheno Phenotypes/pack_years.csvphen \
--analysis-type preprocess \
--fixed_effects Covariates/basic_pack_years.csv \
--fixedEffectNumber 6 \
--thread 24 \
--thread-spawned 24 \
--marker-cache --seed 1 

../../../BayesRRcmd/src/brr \
--data-file Inputs/Methylation/pack_years_17833.csv \
--pheno Phenotypes/pack_years.csvphen \
--fixed_effects Covariates/basic_pack_years.csv \
--fixedEffectNumber 6 \
--analysis-type ppbayes \
--chain-length 25000 \
--burn-in 5000 \
--thin 5 \
--S "0.01,0.1,1.0" \
--mcmc-samples runs/20k_chain_covars/pack_years.csv \
--thread 12 \
--thread-spawned 12 \
--marker-cache --seed 1 > runs/20k_chain_covars/pack_years_20k_chain.out
