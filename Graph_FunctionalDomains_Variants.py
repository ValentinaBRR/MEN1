#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:41:02 2023

@author: valentinaburrai
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
path_graphs = '/Users/valentinaburrai/Data/MEN1/Output/Graphs/'
f_clinvar = 'cv_full.xlsx'

lt_fancd2 = [((219-5), 177)]
lt_disordered = [(465-4, 93)]
lt_basic_acid = [((465-5), 38)]

df = pd.read_excel(path_project+f_clinvar)


#%% 
'''
based on Thakker 2014 (PMID: 23933118)
two variants in his top 9 list could not be matched to ClinVar
'''
#this was the index of the df used to build a univocal list [51,66, 83?, 147, 198, 314, 359, 396, 397]

l_var_id_most_freq = [16693, 265234, 201019, 200997, 200981, 403802, 16692, 200999, 279852]
l_location_hot_spots = df[(df.var_id.isin(l_var_id_most_freq))
                          &
                          (df.aa_location.notnull())].aa_location.tolist()

#%%

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

se = df_ex.groupby('var_type')['var_id'].size()
df_evars = pd.concat([se, se / se.sum() ], axis=1, keys=('count','share'))
#df_vars.loc['Total'] = df_vars.sum().astype(int)
df_evars.sort_values(by='share', ascending=False, inplace=True)
df_evars['tally'] = df_evars.share.cumsum()

sea = df_ex.groupby('aa_location')['aa_location'].size()
df_aa = pd.concat([sea, sea / sea.sum() ], axis=1, keys=('count','share'))
df_aa.reset_index(inplace=True)
df_aa.aa_location.astype(int)

df_men = pd.read_excel(path_project + 'Menin.xlsx')

df_af_graph = pd.merge(df_men, #relative frequency of variants per aa
                 df_aa,
                 how='outer',
                 left_on='position',
                 right_on='aa_location')

df_af_graph_hot_spots = df_af_graph[df_af_graph.aa_location.isin(l_location_hot_spots)].copy()
df_af_graph_other = df_af_graph[~df_af_graph.aa_location.isin(l_location_hot_spots)].copy()

cm = 1/2.54  # centimeters in inches

#%%
fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(15*cm, 8*cm),
                       sharex=True)  #rows, columns, width and height
plt.subplots_adjust(hspace=0)
fig.suptitle('Number of variants by position', fontsize=10)

ax.broken_barh(lt_fancd2, (0, 8),
                  facecolors=('cornflowerblue'),
                  alpha=0.2)
ax.broken_barh(lt_disordered, (0, 8),
                  facecolors=('tomato'),
                  alpha=0.2)
ax.broken_barh(lt_basic_acid,
                  (0, 9),
                  facecolors='mediumaquamarine',
                  alpha=0.2)

ax.bar(df_af_graph_other['position'], df_af_graph_other['count'],
          color='black')
ax.bar(df_af_graph_hot_spots['position'], df_af_graph_hot_spots['count'],
          color='red')

ax.spines['left'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.set_xlabel('Position on the amino acid sequence', fontsize=8)

#ax.set_ylim(5, 35)
ax.set_xlim(0, 610)

plt.xticks(
#    rotation=20,
    ha='right',
    fontsize=8)
plt.yticks(
#    rotation=20,
    ha='right',
    fontsize=8)

ax.tick_params(bottom=True,
                left=False)
#plt.locator_params(axis='x', nbins=2)

ax.annotate('FANCD2 interacting region', (250, 7),
            xytext=(0.63, 0.8), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='cornflowerblue',
            horizontalalignment='right', verticalalignment='top')

ax.annotate('Disordered region', (470, 8),
            xytext=(0.95, 0.8), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='tomato',
            horizontalalignment='right', verticalalignment='top')

ax.annotate('Basic and acid residue region', (470,8),
            xytext=(0.95, 0.95), textcoords='axes fraction',
#            arrowprops=dict(facecolor='gray', shrink=0.95),
            fontsize=8,
            color='mediumaquamarine',
            horizontalalignment='right', verticalalignment='top')


plt.tight_layout()
plt.savefig(fname=path_graphs + 'FunctionalDomainsVariants.png', dpi=300)
