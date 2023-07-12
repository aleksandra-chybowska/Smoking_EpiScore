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

path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesRR/results/runs/complete/"
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

print(sexagnos["PIP"])
fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = sexagnos["PIP"], meta = sexagnos, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "PIP", colors = "custom2")
axes.set_title("Pack Years")
axes.axhline(y=0.95, color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir_ii + "manhattan_py.pdf", dpi=300)
plt.close(fig)