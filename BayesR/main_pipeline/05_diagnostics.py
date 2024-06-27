import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import colorsys
import datatable as dt
import math
from scipy import stats
import statsmodels.stats.multitest as smm
from statsmodels.graphics.tsaplots import plot_pacf
from statsmodels.graphics.tsaplots import plot_acf
import pymc3 as pm3
import arviz as az
sys.path.append('/Cluster_Filespace/Marioni_Group/Ola/Code/general/projects/smoking/BayesR/Elena/plots')      
from plot_funcs import *

# Color palette
palette = "custom2"

def flatten(xss):
    return [x for xs in xss for x in xs]

path = "//Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/results/runs/complete/"
trait = "pack_years_17865_complete"

# output_dir = path + "plots/sex/"
output_dir_ii = path + "plots/explore/"
# output_dir_iii = path + "plots/gene_enrichment/"
# output_dir_iv = path + "plots/masking/"

sexagnos = pd.read_table(path + "summary/" + trait + "_meanbeta_pip.tsv", index_col = 1) 

sigma_sexagnos_file = path + "sigma/" + trait + ".csv"
sigma_sexagnos = pd.read_table(sigma_sexagnos_file, sep = ",")

sigma_sexagnos["varexp"] = sigma_sexagnos.iloc[:,1]/(sigma_sexagnos.iloc[:,0]+sigma_sexagnos.iloc[:,1])

varexplained_sexagnos = pd.read_table(path + "summary/" + trait + "_varianceexplained.tsv", index_col = 0)

prop_varexplained_sexagnos = pd.read_table(path + "summary/" + trait + "_varianceexplained_periteration.tsv")

# Fix chromosome column
sexagnos["chr"] = [int(i.replace("chr", "")) for i in sexagnos["chr"]]


# Convergence
###########################################################

## Parameter values over iterations
i = 0
i_dic = {0 : "sexagnos"}
i_dic_ii = {0 : ""}
i_dic_iii = {0 : "Sex-agnostic"}

# Convergence for each seed (sum of sigmas)
fig, axes = plt.subplots(1, 1, figsize = (7, 3), sharex = True)
col = set_colors(4, palette)
#axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
seed_df["ite"] = range(1,1001)
axes.plot(seed_df["ite"], seed_df["sum"], linestyle = "solid", c = col[i], zorder = 2, label = "20k Iterations", linewidth = 1)

sns.despine(offset=10, trim=True);
axes.set_xlabel("Iteration")
axes.set_ylabel("Sigma Sum")
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain 1") for i in range(0, 1)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.title("Sigma Convergence")
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s.pdf" % i_dic[i], dpi = 300)
plt.close()


# Convergence (proportion variance explained by G)
fig, axes = plt.subplots(1, 1, figsize = (7, 3), sharex = True)
col = set_colors(4, palette)
#axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)

f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
seed_df["ite"] = range(1,1001)
axes.plot(seed_df["ite"], seed_df["prop"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 1)

sns.despine(offset=10, trim=True);
axes.set_xlabel("Iteration")
axes.set_ylabel("Sigma G Prop")
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain 1") for i in range(0, 1)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.title("Sigma G/(Sigma G + Sigma E) Convergence%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_proportion.pdf" % i_dic[i], dpi = 300)
plt.close()

# Convergence (individual sigmas)
fig, axes = plt.subplots(2, 1, figsize = (7, 5), sharex = True)
col = set_colors(4, palette)
#axes[0].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
#axes[1].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)

f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["ite"] = range(1,1001)
axes[0].plot(seed_df["ite"], seed_df["sigmaE"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 1)
axes[1].plot(seed_df["ite"], seed_df["sigmaG[1]"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 1)

sns.despine(offset=10, trim=True);
axes[1].set_xlabel("Iteration")
axes[0].set_ylabel("Sigma E")
axes[1].set_ylabel("Sigma G")
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain 1") for i in range(0, 4)]
# fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.suptitle("Sigma Convergence%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_perparameter.pdf" % i_dic[i], dpi = 300)
plt.close()

# Convergence rolling median plot (sum of sigmas)
fig, axes = plt.subplots(1, 1, figsize = (7, 3))
col = set_colors(4, palette)
#axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
seed_df["ite"] = range(1,1001)
for j in range(0, len(seed_df.index)):
    if j == 0:
        subset = seed_df.loc[seed_df.index[0]]
    else:
        subset = seed_df.loc[seed_df.index[0:j]]
    seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median"] = median(subset["sum"])
axes.plot(seed_df["ite"], seed_df["rolling_median"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 2)

sns.despine(offset=10, trim=True);
axes.set_xlabel("Iteration")
axes.set_ylabel("Median Sigma G + Sigma E")
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain 1") for i in range(0, 4)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.title("Rolling Median Sigma G + Sigma E%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_rollingmedian.pdf" % i_dic[i], dpi = 300)
plt.close()

# Convergence rolling median plot (proportion)
fig, axes = plt.subplots(1, 1, figsize = (7, 3))
col = set_colors(4, palette)
#axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)

f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
seed_df["ite"] = range(1,1001)
for j in range(0, len(seed_df.index)):
    if j == 0:
        subset = seed_df.loc[seed_df.index[0]]
    else:
        subset = seed_df.loc[seed_df.index[0:j]]
    seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median"] = median(subset["prop"])
axes.plot(seed_df["ite"], seed_df["rolling_median"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 2)

sns.despine(offset=10, trim=True);
axes.set_xlabel("Iteration")
axes.set_ylabel("Median Sigma G Prop")
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain 1") for i in range(0, 1)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.title("Rolling Median Sigma G/(Sigma G + Sigma E)%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_rollingmedian_proportion.pdf" % i_dic[i], dpi = 300)
plt.close()

# Convergence rolling median plot (per sigma)
fig, axes = plt.subplots(2, 1, figsize = (7, 5), sharex = True)
col = set_colors(4, palette)
# axes[0].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
# axes[1].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)

f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["ite"] = range(1, 1001)
for j in range(0, len(seed_df.index)):
    if j == 0:
        subset = seed_df.loc[seed_df.index[0]]
    else:
        subset = seed_df.loc[seed_df.index[0:j]]
    seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median_sigmaG"] = median(subset["sigmaG[1]"])
    seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median_sigmaE"] = median(subset["sigmaE"])
axes[0].plot(seed_df["ite"], seed_df["rolling_median_sigmaE"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 2)
axes[1].plot(seed_df["ite"], seed_df["rolling_median_sigmaG"], linestyle = "solid", c = col[i], zorder = 2, label = "Chain 1", linewidth = 2)

sns.despine(offset=10, trim=True);
axes[1].set_xlabel("Iteration")
axes[0].set_ylabel("Median Sigma E")
axes[1].set_ylabel("Median Sigma G")
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain 1") for i in range(0, 4)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.suptitle("Rolling Median Sigma G and Sigma E%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_rollingmedian_perparameter.pdf" % i_dic[i], dpi = 300)
plt.close()

## Histogram per parameter
i = 0
fig, axes = plt.subplots(1, 2, figsize = (6, 3), sharex = True, sharey = True)
col = set_colors(4, palette)
seed = 1
f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
#seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
#seed_df["ite"] = range(1,1001)
sns.histplot(seed_df["sigmaE"], ax = axes[0], color = col[0], label = "Chain %s" % seed, alpha = None, linewidth = 0, edgecolor = "white")
sns.histplot(seed_df["sigmaG[1]"], ax = axes[1], color = col[0], label = "Chain %s" % seed, alpha = None, linewidth = 0, edgecolor = "white")
axes[0].set_xlabel("Sigma E")
axes[1].set_xlabel("Sigma G")
axes[0].set_title("Sigma E")
axes[1].set_title("Sigma G")
axes[1].set_ylabel("Count")
axes[0].set_ylabel("Count")
#axes[0].set_xlim([0,1])
#axes[1].set_xlim([0,1])
sns.despine(offset=10, trim=True);
#plt.suptitle("Sigma Autocorrelation")
plt.tight_layout()
fig.savefig(output_dir_ii + "convergence_smoking_%s_hist.pdf" % i_dic[i], dpi = 300)
plt.close()


## Autocorrelation
fig, axes = plt.subplots(1, 1, figsize = (10, 3), sharex = True, sharey = True)
col = set_colors(4, palette)

seed = 1
f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
seed_df["ite"] = range(1,1001)
plot_acf(seed_df["sum"], ax = axes, color = col[i], label = "Chain %s" % seed, lags = 50, alpha = None)
axes.set_title("Chain %s" % seed)
axes.set_ylabel("Autocorrelation")
axes.set_xlabel("Lag")
axes.set_ylim([-0.1,1.1])

sns.despine(offset=10, trim=True);
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.suptitle("Sigma Sum Autocorrelation%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_autocorrelation.pdf" % i_dic[i], dpi = 300)
plt.close()

# Convergence for each seed (proportion)
fig, axes = plt.subplots(1, 1, figsize = (10, 3), sharex = True, sharey = True)
col = set_colors(4, palette)
seed = 1
f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
seed_df["ite"] = range(1,1001)
plot_acf(seed_df["prop"], ax = axes, color = col[i], label = "Chain %s" % seed, lags = 50, alpha = None)
axes.set_title("Chain %s" % seed)
axes.set_ylabel("Autocorrelation")
axes.set_xlabel("Lag")
axes.set_ylim([-0.1,1.1])

sns.despine(offset=10, trim=True);
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.suptitle("Sigma G/(Sigma G + Sigma E) Autocorrelation%s" % i_dic_ii[i])
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_autocorrelation_proportion.pdf" % i_dic[i], dpi = 300)
plt.close()


# Convergence for each seed (individual sigmas)
fig, axes = plt.subplots(2, 1, figsize = (10, 5), sharex = True, sharey = True)
col = set_colors(4, palette)
i = 0
seed = 1
f = sigma_sexagnos_file
seed_df = pd.read_table(f, sep = ",")
seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
seed_df["ite"] = range(1,1001)
plot_acf(seed_df["sigmaE"], ax = axes[0], color = col[i], label = "Chain %s" % seed, lags = 50, alpha = None)
plot_acf(seed_df["sigmaG[1]"], ax = axes[1], color = col[i], label = "Chain %s" % seed, lags = 50, alpha = None)
#axes[0,s].acorr(seed_df["sigmaE"], color = col[s], label = "Chain %s" % seed, maxlags = 50)
#axes[1,s].acorr(seed_df["sigmaG[1]"], color = col[s], label = "Chain %s" % seed, maxlags = 50)
axes[0].set_title("Chain %s" % seed)
axes[0].set(xlabel = None)
#axes[1,s].set_xlabel("Lag")
axes[1].set(title = None)

axes[1].set_ylabel("Sigma G")
axes[0].set_ylabel("Sigma E")
sns.despine(offset=10, trim=True);
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
#fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.suptitle("Sigma Autocorrelation%s" % i_dic_ii[i])
fig.supylabel("Autocorrelation")
fig.supxlabel("Lag")
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_smoking_%s_autocorrelation_perparameter.pdf" % i_dic[i], dpi = 300)
plt.close()

## Geweke diagnostics
# i = 0
# fig, axes = plt.subplots(3, 1, figsize = (6, 5))
# col = set_colors(4, palette)

# seed = 1
# f = sigma_sexagnos_file
# seed_df = pd.read_table(f, sep = ",")
# geweke_sigmaG = pm3.geweke(seed_df["sigmaG[1]"], intervals = 1, first=0.1, last=0.5)
# geweke_sigmaE = pm3.geweke(seed_df["sigmaE"], intervals = 1, first=0.1, last=0.5)
# print(geweke_sigmaG[0][1])
# print(geweke_sigmaE[0][1])
# axes[i].plot(geweke_sigmaG[0][1], 0, marker = 'o', label = "Chain %s" % seed, color = col[s], zorder = 2)
# axes[i].plot(geweke_sigmaE[0][1], 1, marker = 'o', label = "Chain %s" % seed, color = col[s], zorder = 2)
# axes[i].set_xlim([-5, 5])
# axes[i].set_title(i_dic_iii[i])
# axes[i].set_yticks([0, 1])
# axes[i].set_yticklabels(["Sigma G", "Sigma E"])
# axes[i].set_ylim([-0.1, 1.1])
# axes[i].axvline(-2, color = 'grey', alpha = 0.7, linestyle = "dashed", linewidth = 1.5, zorder = 1)
# axes[i].axvline(2, color = 'grey', alpha = 0.7, linestyle = "dashed", linewidth = 1.5, zorder = 1)

# axes[2].set_xlabel("Z")
# fig.suptitle("Geweke Z Diagnostic")
# plt.tight_layout()
# sns.despine(offset=10, trim=True);
# elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
# # fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
# fig.subplots_adjust(right=0.87) 
# fig.savefig(output_dir_ii + "convergence_smoking_geweke_perparameter.pdf", dpi = 300)
# plt.close()

# ## Effective sample size
# i = 0
# ess_df = pd.DataFrame(index = flatten(list(i_dic[i])), columns = ["Parameter"] + ["Chain 1"])
# ess_df["Parameter"] = ["Sigma G", "Sigma E", "Sigma Sum", "Sigma Prop"]

# seed = 1
# f = sigma_sexagnos_file
# seed_df = pd.read_table(f, sep = ",")
# seed_df["sum"] = seed_df["sigmaE"] + seed_df["sigmaG[1]"]
# seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
# ess_G = az.ess(np.array(seed_df["sigmaG[1]"]))
# ess_E = az.ess(np.array(seed_df["sigmaE"]))
# ess_sum = az.ess(np.array(seed_df["sum"]))
# ess_prop = az.ess(np.array(seed_df["prop"]))
# ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma G"), "Chain %s" % seed] = ess_G
# ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma E"), "Chain %s" % seed] = ess_E
# ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma Sum"), "Chain %s" % seed] = ess_sum
# ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma Prop"), "Chain %s" % seed] = ess_prop

# ess_df.to_csv(output_dir_ii + "ess_acrosschains.tsv", sep = "\t", index_label = "Model", na_rep = "NA")

#manhattan

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = sexagnos["PIP"], meta = sexagnos, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "PIP", colors = "custom2")
axes.set_title("Alcohol Consumption")
axes.axhline(y=0.95, color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir_ii + "manhattan_alcoholconsumption_sexagnos.pdf", dpi=300)
plt.close(fig)