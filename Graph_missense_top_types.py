#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:56:44 2023

@author: valentinaburrai
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
path_graphs = '/Users/valentinaburrai/Data/MEN1/Output/Graphs/'
f_cv = 'cv_full.xlsx'
#f_ensembl = 'Ensembl_exons_2023_05_01_ENST00000450708.csv'
'''
loads the file extracted from Clinvar into a pandas dataframe
'''


df = pd.read_excel(path_project + f_cv)
l_introns = df[
    (df.var_type.notnull())
    &
    ((df.var_type.str.contains('intron'))
     |
     (df.var_type.str.contains('splic')))
    ].var_id.tolist()
l_nulls = df[
    (df.var_type.isnull())
    ].var_id.tolist()
l_exons = df[
    (~df.var_id.isin(l_introns))
    &
    (~df.var_id.isin(l_nulls))
             ].var_id.tolist()
df_ex = df[df.var_id.isin(l_exons)].copy()

df_m = df[df.missense_type.notnull()].copy()
s1 = df_m.missense_type.value_counts()
l_tops = s1.index[s1>=5].tolist()

df_m['missense_top_types'] = np.where(df_m.missense_type.isin(l_tops),
                                      df_m.missense_type,
                                      'other')
l_tops.append('other')


l_subs = df_ex.missense_type[df_ex.missense_type.notnull()].unique().tolist()

l_colours = sns.color_palette('colorblind',len(l_tops))

d_markers = {
    'neutral_to_neutral' : '<',
    'neutral_to_charged_positive': '>',
    'polar_uncharged_to_aromatic': 'd',
    'neutral_to_charged_negative': 's',
    'aromatic_to_polar_uncharged': 'p',
    'charged_negative_to_charged_positive': 'X',
    'neutral_to_polar_uncharged': '^',
    'other': 'v'
    }



'''
generates a standardised figure and graph
for use in word documents'''
cm = 1/2.54  # centimeters in inches
fig, axs = plt.subplots(nrows=len(l_tops),
                       ncols=1,
                       figsize=(15*cm, 8*cm),
                       sharex=True)  #rows, columns, width and height
plt.subplots_adjust(hspace=0)
fig.suptitle('Most frequent types of missense substitutions', fontsize=10)
for i, ax in enumerate(axs.flat):

    sns.stripplot(
        data=df_m[['missense_top_types', 'aa_location']][df_m.missense_top_types == l_tops[i]],
        ax=ax,
        x="aa_location",
        y="missense_top_types",
#        hue="var_type",
#        kind="swarm",
        size=2,
#        edgecolor=l_colours[i],
        color=l_colours[i],
        marker = d_markers[l_tops[i]]
                    )
    ax.set_title(str(l_tops[i]).replace('_',' ').capitalize(),
                 color=l_colours[i],
                 fontsize=8,
                 y=0.7, x=0.8)
    ax.grid(visible=True, which='major', axis='x')
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlim(0, 610)
    plt.xticks(
    #    rotation=20,
        ha='right',
        fontsize=8)

    ax.set_xlabel('Position on the amino acid sequence', fontsize=8)
    ax.set_ylabel('')
    ax.tick_params(bottom=False,
                    left=False)


plt.tight_layout()
plt.savefig(fname=path_graphs + 'Missense.png', dpi=300)
