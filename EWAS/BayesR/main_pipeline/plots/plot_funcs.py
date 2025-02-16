#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import colorsys
import datatable as dt
import math

# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

font = "Galvji"
path = '/Cluster_Filespace/Marioni_Group/Elena/other/fonts/Galvji/Galvji-01.ttf'
fe = matplotlib.font_manager.FontEntry(
    fname= path,
    name=font)
matplotlib.font_manager.fontManager.ttflist.insert(0, fe)
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = font

def lighten_color(color, amount=1.2):  
    # --------------------- SOURCE: @IanHincks ---------------------
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def set_colors(N, name = "default"):
    colors = []
    if name == "default":
        colors = sns.color_palette("Set2", 6)
    elif name == "custom1":
        colors = ["#DD1155", "#7D82B8", "#FCBA04", "#63C7B2", "#45364B"] #Hotpink, lilacky-blue, teal, yellow, dark purple
    elif name == "custom2":
        colors = ["#47C3F4", "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B", "07B1D8"] #Lightblue, orange, teal, yellow, red, dark blue, midtone blue
    colors = (colors*int(ceil(N/len(colors))))[:N]
    return colors

def cpg_trajectory(df, cpg, gene_name, output_dir, cpg_col = "CpG", age_col = "Age", color = "#07B1D8", wave = "w3"):
    # Plot - all together
    fig, axes = plt.subplots(1, 2, figsize = (6.5, 3.5), sharey=True, gridspec_kw=dict(width_ratios=[4,1]))
    sns.regplot(ax = axes[0], data = df, x = age_col, y = cpg_col, color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(color)}, lowess = True)
    sns.kdeplot(ax = axes[1], data = df, y = cpg_col, shade = True, color = color, linewidth = 2)
    if "beta" not in wave:
        axes[0].set_ylabel("CpG M-Value")
    else:
        axes[0].set_ylabel("CpG Beta-Value")
        axes[0].set_ylim([0,1])
    axes[0].set_xlabel("Age")
    plt.suptitle("%s | %s" % (cpg, gene_name))
    plt.tight_layout()
    plt.savefig(output_dir + "%s_%s_trajectory_%s.pdf" % (cpg, gene_name, wave), dpi = 300)
    plt.close(fig)

def cpg_trajectory_nodens(df, cpg, gene_name, output_dir, cpg_col = "CpG", age_col = "Age", color = "#07B1D8", wave = "w3"):
    # Plot - all together
    fig, axes = plt.subplots(1, 1, figsize = (5, 3.5))
    sns.regplot(ax = axes, data = df, x = age_col, y = cpg_col, color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':"black"}, lowess = True)
    sns.kdeplot(data = df, x = age_col, y = cpg_col, color = lighten_color(color), levels = 4, thresh=.1, linewidth=0.5)
    if "beta" not in wave:
        axes.set_ylabel("CpG M-Value")
    else:
        axes.set_ylabel("CpG Beta-Value")
        axes.set_ylim([0,1])
    axes.set_xlabel("Age")
    axes.set_title("%s | %s" % (cpg, gene_name))
    plt.tight_layout()
    plt.savefig(output_dir + "%s_%s_trajectory_%s_nodens.pdf" % (cpg, gene_name, wave), dpi = 300)
    plt.close(fig)

def genomic_pos(pos, col_chr, col_bp):
    l = pos.groupby(col_chr)[col_bp].max()
    offset = l.cumsum() - l
    return pos[col_bp] + list(offset.loc[pos[col_chr]])

def chromosome_colors(N, name='custom2'):
    colors = []
    if name == "default":
        colors = sns.color_palette("Set2", 6)
    elif name == "custom1":
        colors = ["#DD1155", "#7D82B8", "#FCBA04", "#63C7B2", "#45364B"] #Hotpink, lilacky-blue, teal, yellow, dark purple
    elif name == "custom2":
        colors = ["#47C3F4", "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B"] #Lightblue, orange, teal, yellow, red, dark blue, midtone blue, very light blue, 
    colors = (colors*int(math.ceil(N/len(colors))))[:N]
    return colors

def manhattan(values,meta,col_chr=None,col_bp=None,ylabel=None,ax=None,ylim=None,colors=None,annotate=None):
    meta = meta.loc[values.index].copy()
    meta['genomic_position'] = genomic_pos(meta,col_chr,col_bp)
    data = pd.concat([meta,values.to_frame('value')],axis=1)
    Nchr = len(data[col_chr].unique())

    if colors is None:
        colors = chromosome_colors(Nchr)
    elif type(colors) is str:
        colors = chromosome_colors(Nchr, name=colors)
    
    ticks = {}
    for c,g in data.groupby(col_chr):
        bounds = [g.genomic_position.min(),g.genomic_position.max()]
        ticks[c] = bounds[0] + (bounds[1]-bounds[0])/2
    print(ticks)
    
    if ax is None:
        fig, ax = subplots(figsize=(15,8))
        
    which = data.index if ylim is None else (ylim[0] <= data.value) & (data.value <= ylim[1])
    i = 0
    for c,g in data.loc[which].groupby(col_chr):
        if c != "X":
            i += 1
            color = colors[i-1]
        else:
            color = "slategrey"
        print(g.genomic_position)
        print(g.value)
        sns.scatterplot(data=g, x="genomic_position", y="value", color=color, edgecolor=None, marker ='.', rasterized = True)
    
    ax.set_xticks(list(ticks.values()))
    ax.set_xticklabels(list(ticks.keys()))
    ax.tick_params()
    ax.set_xlabel('Genomic Position')
    if not ylabel is None:
        ax.set_ylabel(ylabel)

def densityplot(df, varname, ax = None, font = "Galvji", col = "custom1"):
    df = df[~df[varname].isnull()].copy()
    palette = set_colors(1, col)
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font
    sns.kdeplot(ax = ax, data = df, x = varname, shade = True, color = palette, linewidth = 2)
    ax.set_xlabel(varname.capitalize())
    ax.set_ylabel('Density')

def densityplot_batches(df, varname, ax = None, batch_column = "batch", batches = [1,3,4], font = "Galvji", col = "custom1"):
    df = df[~df[varname].isnull()].copy()
    palette = set_colors(len(batches), col)
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font
    i = 0
    for b in batches:
        df_b = df[df[batch_column] == b]
        sns.kdeplot(ax = ax, data = df_b, x = varname, shade = True, label = "Wave %s" % b, color = palette[i], linewidth = 2)
        ax.set_xlabel(varname.capitalize())
        ax.set_ylabel('Density')
        i += 1

def boxplot(df, varname, ax = None, font = "Galvji", col = "custom1"):
    df = df[~df[varname].isnull()].copy()
    df["ph"] = "" #Placeholder
    palette = set_colors(1, col)
    sns.set_palette(sns.color_palette(palette))
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font
    sns.boxplot(ax = ax, data = df, x = "ph", y = varname, linewidth = 1, saturation = 1)  
    ax.set_xlabel("All waves")
    ax.set_ylabel(varname.capitalize())

    for i,artist in enumerate(ax.artists):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = lighten_color(artist.get_facecolor(), 1.2)
        artist.set_edgecolor(col)    

        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        for j in range(i*6,i*6+6):
            line = ax.lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)

def boxplot_batches(df, varname, ax = None, batch_column = "batch", batches = [1,3,4], font = "Galvji", col = "custom1"):
    df = df[~df[batch_column].isnull()].copy()
    df = df[~df[varname].isnull()].copy()
    df.loc[:,batch_column] = df.loc[:,batch_column].astype(int).copy()
    df = df.sort_values(by = batch_column, ascending = True)
    
    palette = set_colors(len(batches), col)
    sns.set_palette(sns.color_palette(palette))
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font
    sns.boxplot(ax = ax, data = df, x = batch_column, y = varname, linewidth = 1, saturation = 1)  
    ax.set_xlabel("Wave")
    ax.set_ylabel(varname.capitalize())

    for i,artist in enumerate(ax.artists):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = lighten_color(artist.get_facecolor(), 1.2)
        artist.set_edgecolor(col)    

        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        for j in range(i*6,i*6+6):
            line = ax.lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)

def histogram(df, varname, ax = None, font = "Galvji", col = "custom1", nr_bins = 20):
    df = df[~df[varname].isnull()].copy()
    palette = set_colors(1, col)
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font
    
    sns.histplot(ax = ax, data = df, x = varname, color = palette, alpha=0.6, edgecolor=palette, linewidth = 0, bins = nr_bins)
    ax.set_xlabel(varname.capitalize())
    ax.set_ylabel('Frequency')

def histogram_batches(df, varname, ax = None, batch_column = "batch", batches = [1,3,4], font = "Galvji", col = "custom1"):
    df = df[~df[varname].isnull()].copy()
    palette = set_colors(len(batches), col)
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font
    
    i = 0
    for b in batches:
        df_b = df[df[batch_column] == b]
        sns.histplot(ax = ax, data = df_b, x = varname, label = "Wave %s" % b, color = palette[i], alpha=0.6, edgecolor=palette[i], linewidth = 0)
        ax.set_xlabel(varname.capitalize())
        ax.set_ylabel('Frequency')
        i += 1

def countplot(df, varname, ax = None, font = "Galvji", col = "custom1", legend = False):
    df = df[~df[varname].isnull()].copy()
    palette = set_colors(1, col)
    sns.set_style("white", {'legend.frameon':True})
    sns.set_palette(sns.color_palette(palette))
    matplotlib.rcParams['font.family'] = font

    sns.countplot(ax = ax, data = df, x = varname, linewidth = 0, saturation = 1, color = palette[0])
    #for patch in ax.patches:
    #    col = lighten_color(patch.get_facecolor(), 1.2)
    #    patch.set_edgecolor(col)
   
    ax.set_xlabel(varname.capitalize())
    ax.set_ylabel("Count")

    if legend == False:
        try:
            ax.get_legend().remove()
        except:
            print("\t\tNote: No legend to remove.")

def countplot_batches(df, varname, ax = None, batch_column = "batch", batches = [1,3,4], font = "Galvji", col = "custom1", legend = True):
    df = df[~df[batch_column].isnull()].copy()
    df.loc[:,batch_column] = df.loc[:,batch_column].astype(int).copy()
    df = df.sort_values(by = batch_column, ascending = True)
    df.loc[:,batch_column] = df.loc[:,batch_column].apply(lambda x: "Wave %s" % x).copy()
    
    palette = set_colors(len(batches), col)
    sns.set_style("white", {'legend.frameon':True})
    sns.set_palette(sns.color_palette(palette))
    matplotlib.rcParams['font.family'] = font

    sns.countplot(ax = ax, data = df, x = varname, hue = batch_column, hue_order = list(map(lambda x: "Wave %s" % x, batches)), linewidth = 0, saturation = 1)
    #for patch in ax.patches:
    #    col = lighten_color(patch.get_facecolor(), 1.2)
    #    patch.set_edgecolor(col)
   
    ax.set_xlabel(varname.capitalize())
    ax.set_ylabel("Count")

    if legend == False:
        try:
            ax.get_legend().remove()
        except:
            print("\t\tNote: No legend to remove.")

def proportionplot(df, varname, ax = None, font = "Galvji", col = "custom1", legend = True):
    palette = set_colors(1, col)
    sns.set_style("white", {'legend.frameon':True})
    matplotlib.rcParams['font.family'] = font

    df = df[~df[varname].isnull()]
    g_df = pd.DataFrame(df[varname].value_counts(normalize=True))
    g_df.columns = ["Percentage"]
    g_df[varname] = g_df.index.values.tolist().copy()
    
    sns.barplot(ax = ax, x=varname, y='Percentage', data=g_df, saturation = 1, color = palette[0], linewidth = 0)
    #for patch in ax.patches:
    #    col = lighten_color(patch.get_facecolor(), 1.2)
    #    patch.set_edgecolor(col)
   
    ax.set_xlabel(varname.capitalize())
    ax.set_ylabel("Percentage")

    if legend == False:
        try:
            ax.get_legend().remove()
        except:
            print("\t\tNote: No legend to remove.")

def proportionplot_batches(df, varname, ax = None, batch_column = "batch", batches = [1,3,4], font = "Galvji", col = "custom1", legend = True):
    palette = set_colors(len(batches), col)
    sns.set_style("white", {'legend.frameon':True})
    sns.set_palette(sns.color_palette(palette))
    matplotlib.rcParams['font.family'] = font

    df = df[~df[batch_column].isnull()]
    df = df[~df[varname].isnull()]
    df = df.sort_values(by = batch_column, ascending = True)
    df.loc[:,batch_column] = df.loc[:,batch_column].astype(int).copy()
    df.loc[:,batch_column] = df.loc[:,batch_column].apply(lambda x: "Wave %s" % x).copy()
    y,x = varname,batch_column
    g_df = df.groupby(x)[y].value_counts(normalize=True)
    g_df = g_df.rename('Percentage').reset_index()
    
    sns.barplot(ax = ax, x=varname, y='Percentage', hue=batch_column, hue_order = list(map(lambda x: "Wave %s" % x, batches)), data=g_df, saturation = 1, linewidth = 0)
    #for patch in ax.patches:
    #    col = lighten_color(patch.get_facecolor(), 1.2)
    #    patch.set_edgecolor(col)
   
    ax.set_xlabel(varname.capitalize())
    ax.set_ylabel("Percentage")

    if legend == False:
        try:
            ax.get_legend().remove()
        except:
            print("\t\tNote: No legend to remove.")

