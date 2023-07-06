#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 21:45:19 2023

@author: valentinaburrai

This codes two heatmaps for substitutions of amino acids in missense
 pathological mutations and for substitutions of amono acid by type in
 missense pathological mutations.

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
path_graphs = '/Users/valentinaburrai/Data/MEN1/Output/Graphs/'
f_clinvar = 'cv_full.xlsx'


df = pd.read_excel(path_project+f_clinvar)

df = df[df.var_class == 'germline'].copy()

df_m = df[df.var_type == 'missense'].copy()

shares_aa = pd.DataFrame(df_m.aa_wild_type.value_counts())
shares_aa['pct_share'] = shares_aa.aa_wild_type/shares_aa.aa_wild_type.sum()
shares_aa['cumulative_share'] = shares_aa.pct_share.cumsum()

shares_aa2 = pd.DataFrame(df_m.aa_mutated.value_counts())
shares_aa2['pct_share'] = shares_aa2.aa_mutated/shares_aa2.aa_mutated.sum()
shares_aa2['cumulative_share'] = shares_aa2.pct_share.cumsum()

shares_at = pd.DataFrame(df_m.aa_wt_type.value_counts())
shares_at['pct_share'] = shares_at.aa_wt_type/shares_at.aa_wt_type.sum()
shares_at['cumulative_share'] = shares_at.pct_share.cumsum()

shares_at2 = pd.DataFrame(df_m.aa_m_type.value_counts())
shares_at2['pct_share'] = shares_at2.aa_m_type/shares_at2.aa_m_type.sum()
shares_at2['cumulative_share'] = shares_at2.pct_share.cumsum()


df_graph1 = pd.crosstab(
    index=df_m.aa_wild_type,
#    margins=True,
    columns=df_m.aa_mutated
    )

df_graph2 = pd.crosstab(
    index=df_m.aa_wt_type,
    margins=True,
    normalize=True,
    columns=df_m.aa_m_type
    )
l_labels = df_graph2.columns.tolist()
cm = 1/2.54  # centimeters in inches
#%%
'''
heatmap of AA pairwise substitutions.

Unfortunately, this release of matplolib
(3.7.0) has an issue with heatmaps and does not allow control of the font size of the 
tick labels. I left it in there as it may get fixed at the next matplotlib release 
'''
fig1, ax1 = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(14*cm, 14*cm),
                       )  #rows, columns, width and height
plt.subplots_adjust(hspace=0)
fig1.suptitle('Missense amino acid substitutions', fontsize=12)

sns.heatmap(
    df_graph1,
    cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True),
    cbar_kws={'format': '%d'},
    ax=ax1,
#    annot_kws={"fontsize":6},
#    mask=True
    )
ax1.set_xlabel('Mutated amino acid', fontsize=9)
#ax1.set_xticklabels(
#    [str(x).replace('_',' ').capitalize() for x in df_graph1.columns.tolist()],
#    fontsize=7
#                    )
#ax1.set_xticks(None)

ax1.set_ylabel('Wild-type amino acid', fontsize=9)
#ax1.set_yticklabels(
#    [str(x).replace('_',' ').capitalize() for x in df_graph1.columns.tolist()],
#    fontsize=7
#                    )
#ax1.set_yticks(None)

plt.tick_params(bottom=False,
                left=False)

plt.tight_layout()
plt.savefig(fname=path_graphs + 'MissenseAAHeatmap.png', dpi=300)

#%%
'''
heatmap of types of amino acid pairwise substitutions
'''
fig2, ax2 = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(14*cm, 14*cm),
                       )  #rows, columns, width and height

plt.subplots_adjust(hspace=0)
fig2.suptitle('Missense substitutions by amino acid type', fontsize=10)

sns.heatmap(
    df_graph2,
    cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True),
    cbar_kws={'format': '%d'},
    ax=ax2,
#    mask=True
    )

ax2.set_xlabel('Mutated amino acid type', fontsize=9)
ax2.set_xticklabels(
    [str(x).replace('_',' ').capitalize() for x in df_graph2.columns.tolist()],
    fontsize=8
                    )
ax2.set_ylabel('Wild-type amino acid type', fontsize=9)
ax2.set_yticklabels(
    [str(x).replace('_',' ').capitalize() for x in df_graph2.columns.tolist()],
    fontsize=7
                    )


plt.tick_params(bottom=False,
                left=False)

plt.tight_layout()
plt.savefig(fname=path_graphs + 'MissenseTypeHeatmap.png', dpi=300)
