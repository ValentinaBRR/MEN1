#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 16:43:50 2023

@author: valentinaburrai
"""

import pandas as pd
import matplotlib.pyplot as plt

'''
sets names for directories and file
'''
path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
path_graphs = '/Users/valentinaburrai/Data/MEN1/Output/Graphs/'

f_clinvar = 'cv_full.xlsx'

'''
loads Clinvar data tagged with Ensembl, UniProt and NextProt information
and tagged for type of variant
'''
df = pd.read_excel(path_project + f_clinvar)
df = df[df.var_class == 'germline']
df.condition.value_counts()

df_c = df[['var_id', 'condition']].join(df.condition.str.split('|', expand=True))
df_c.drop('condition', axis='columns', inplace=True)
df_c = df_c.melt(id_vars='var_id')
df_c.drop('variable', axis='columns', inplace=True)
df_c.dropna(inplace=True)

df_c = df_c.sort_values(by=['var_id', 'value'], ascending=[True, True])
df_c.value.value_counts()
l_generic = df_c.var_id[df_c.value.isin(['Multiple endocrine neoplasia, type 1',
                                         'Hereditary cancer-predisposing syndrome'])
                        ].unique().tolist()
l_unknown = df_c.var_id[df_c.value.isin(['not provided',
                                         'not specified'])
                        ].unique().tolist()
l_specific = df_c.var_id[(~df_c.value.isin(['not provided',
                                         'not specified',
                                         'Multiple endocrine neoplasia, type 1',
                                         'Hereditary cancer-predisposing syndrome']))
                        ].unique().tolist()

df_s = df[df.var_id.isin(l_specific)].copy()
df_s.var_starts = -abs(df_s.var_starts)


df_graph = pd.merge(df_s,
                    df_c,
                    how='left',
                    on='var_id')
df_graph = df_graph[df_graph.value != 'Multiple endocrine neoplasia, type 1'].copy()
l_conditions = df_graph.value.unique().tolist()
#%%
import numpy as np

n=20
colour = iter(plt.cm.rainbow(np.linspace(0, 1, n)))
l_colour = ['darkgrey',
 'tan',
 'powderblue',
 'lightcoral',
 'forestgreen',
 'lightskyblue',
 'plum',
 'pink',
 'olivedrab',
 'tomato',
 'mediumseagreen',
 'orange',
 'royalblue',
 'rosybrown',
 'orchid',
 'mistyrose',
 'lightsteelblue',
 'mediumpurple',
'gold',
'chartreuse'
 ]
'''
for i in range(n):
    print(l_colour[i])
   c = next(color)
   plt.plot(x, y, c=c)
'''
cm = 1/2.54
fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(15*cm, 8*cm),
                       sharex=True)  #rows, columns, width and height
plt.subplots_adjust(hspace=0)
fig.suptitle('Variants with detailed conditions', fontsize=10)

for v in l_conditions:
    ax.stem(df_graph[df_graph.value == v].var_starts,
            df_graph[df_graph.value == v].value,
            linefmt='grey',  #colour five of the series, solid line
            markerfmt=l_colour[l_conditions.index(v)],  #colour six of the series, circle
            basefmt=" ",  #no baseline
            bottom=-1)

'''

ax.stem(df_graph[df_graph.value == 'Angiofibroma, somatic'].var_starts,
        df_graph[df_graph.value == 'Angiofibroma, somatic'].value,
        linefmt='grey',  #colour five of the series, solid line
        markerfmt='C6o',  #colour six of the series, circle
        basefmt=" ",  #no baseline
        bottom=-2)

ax.stem(df_graph[df_graph.value == 'Somatotroph adenoma'].var_starts,
        df_graph[df_graph.value == 'Somatotroph adenoma'].value,
        linefmt='grey',  #colour five of the series, solid line
        markerfmt='C7o',  #colour six of the series, circle
        basefmt=" ",  #no baseline
        bottom=-2)
ax.stem(df_graph.var_starts,
        df_graph.value,
        linefmt='C5-',  #colour five of the series, solid line
        markerfmt='C6o',  #colour six of the series, circle
        basefmt=" ",  #no baseline
        bottom=-2)
'''
ax.set_xlim(-64810551-400, -64803516+400)
ax.set_xlabel('Chromosome 11q, minus strand', fontsize=9)
ax.set_yticklabels([])
plt.xticks(rotation=20)
#ax.set_xticklabels(
#    range(-64810551-400, -64803516+400),
#    bin(5),
#    rotation=-20, 
#    ha="right",  
#    rotation_mode="anchor")
plt.ticklabel_format(axis='x',
                     
                     style='plain')


ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.tick_params(bottom=False,
               left=False)
ax.annotate('Somatotroph adenoma', (-64807000, 5),
            xytext=(0.82, 0.98), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='orchid',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Adrenocortical adenoma', (-64807000, 5),
            xytext=(0.8, 0.93), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='rosybrown',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Ependymoma', (-64807000, 5),
            xytext=(0.41, 0.87), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='royalblue',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Thyroid adenoma', (-64807000, 5),
            xytext=(0.38, 0.80), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='orange',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Hypertensive disorder', (-64807000, 5),
            xytext=(0.38, 0.74), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='mediumseagreen',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('GI stromal tumour', (-64807000, 5),
            xytext=(0.38, 0.68), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='tomato',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Diabetes mellitus', (-64807000, 5),
            xytext=(0.38, 0.61), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='olivedrab',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Chronic Diarrhea', (-64807000, 5),
            xytext=(0.38, 0.55), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='pink',
            horizontalalignment='right', verticalalignment='top')
'''
ax.annotate('Angiofibroma, somatic', (-64807000, 5),
            xytext=(0.65, 0.52), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='pink',
            horizontalalignment='right', verticalalignment='top')
'''
ax.annotate('Pancreatic insulin-producing NET', (-64807000, 5),
            xytext=(0.455, 0.49), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='plum',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Parathyroid gland adenoma', (-64807000, 5),
            xytext=(0.41, 0.43), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='lightskyblue',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Medullary thyroid carcinoma', (-64807000, 5),
            xytext=(0.42, 0.36), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='forestgreen',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Calcium nephrolithiasis', (-64807000, 5),
            xytext=(0.375, 0.31), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='lightcoral',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Abnormal CCC', (-64807000, 5),
            xytext=(0.295, 0.25), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='powderblue',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Primary hyperparatyroidism', (-64807000, 5),
            xytext=(0.415, 0.185), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='tan',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('Metastatic PNET', (-64807000, 5),
            xytext=(0.35, 0.125), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=7,
            color='darkgrey',
            horizontalalignment='right', verticalalignment='top')



plt.tight_layout()

plt.savefig(fname=path_graphs + 'SpecificConditions.png', dpi=300)



