#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:59:28 2023

@author: valentinaburrai
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
f_clinvar = 'cv_full.xlsx'
f_ensembl = 'Ensembl_exons_2023_05_01_ENST00000450708.csv'
path_graphs = '/Users/valentinaburrai/Data/MEN1/Output/Graphs/'


df_e = pd.read_csv(path_project + f_ensembl, sep=',', header=0)
df_e['Start'] = df_e['Start'].str.replace(',', '')
df_e['End'] = df_e['End'].str.replace(',', '')
df_e['Start'] = df_e['Start'].astype(float).astype('Int64')
df_e['End'] = df_e['End'].astype(float).astype('Int64')
df_e['ensembl_code'] = df_e['Exon / Intron'][df_e['Exon / Intron'].str.startswith('ENSE')]
df_e['exon'] = 'exon_' + df_e['No.']

df_e['Exon / Intron'] = np.where(
    df_e['Exon / Intron'].str.startswith('ENSE'),
    df_e['exon'],
    df_e['Exon / Intron'])
df_e.rename({'Exon / Intron':'exon_intron'}, axis='columns', inplace=True)
df_e.exon_intron = df_e.exon_intron.str.lower(
    ).str.replace(' ', '_')

dfe = df_e.iloc[1:-1,1:-1].reset_index().copy()
dfe.rename(columns={'index': 'order'}, inplace=True)
dfe.drop(['Start Phase', 'End Phase', 'Sequence'], axis='columns', inplace=True)

dfe.columns = dfe.columns.str.replace(' ', '_').str.lower()

dfe.length = dfe.length.str.replace(',', '')
dfe.length = dfe.length.astype(float).astype(int)
dfe['colour'] = ''

l_colours = sns.color_palette('colorblind', 2)
dfe.colour[dfe.exon_intron.str.startswith('e')] = 'cornflowerblue'
dfe.colour[~dfe.exon_intron.str.startswith('e')] = 'tomato'

#%%

cm = 1/2.54  # centimeters in inches


l_colours = sns.color_palette('colorblind', 2)

dfee = dfe.copy()
dfee.start = -abs(dfee.start)
dfee.end = -abs(dfee.end)

graph_starts=dfee.start.min()


lt_exons =list(zip(*map(dfee[dfee.exon_intron.str.startswith('e')].get,
                        ['start', 'length'])))

lt_introns =list(zip(*map(dfee[~dfee.exon_intron.str.startswith('e')].get,
                        ['start', 'length'])))

fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(15*cm, 4*cm),
#                       sharex=True
                       )  #rows, columns, width and height
fig.suptitle('MEN1 transcript NM_001370259.2 ', fontsize=10)

ax.broken_barh(lt_introns, (10, 9),facecolors='tomato')
ax.broken_barh(lt_exons, (20, 9), facecolors='cornflowerblue')

#ax.set_ylim(5, 35)
ax.set_xlim(dfee.start.min()-400, dfee.end.max()+400)
plt.xticks(rotation=20, ha='right', fontsize=8)
#plt.locator_params(axis='x', nbins=2)
ax.set_xlabel('Chromosome 11q, minus strand', fontsize=9)

ax.set_yticks([15, 25], labels=['introns', 'exons'])

ax.grid(False)  # Make grid lines visible
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.annotate('1', (-64810551, 25),
            xytext=(0.06, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')

ax.annotate('2', (-64810551, 25),
            xytext=(0.14, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')

ax.annotate('3', (-64810551, 25),
            xytext=(0.382, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('4', (-64810551, 25),
            xytext=(0.43, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('5', (-64810551, 25),
            xytext=(0.485, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('6', (-64810551, 25),
            xytext=(0.505, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('7', (-64810551, 25),
            xytext=(0.6, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('8', (-64810551, 25),
            xytext=(0.675, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('9', (-64810551, 25),
            xytext=(0.75, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')
ax.annotate('10', (-64810551, 25),
            xytext=(0.88, 0.9), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='black',
            horizontalalignment='right', verticalalignment='top')

plt.ticklabel_format(axis='x', style='plain')

plt.tight_layout()
plt.savefig(fname=path_graphs + 'Transcript.png', dpi=300)

