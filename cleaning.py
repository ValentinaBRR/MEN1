#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:55:57 2023

@author: valentinaburrai
"""
import numpy as np
import pandas as pd


'''
loads the file extracted from Clinvar into a pandas dataframe
'''
path_project = '/Users/valentinaburrai/Downloads/'
f_risultati = 'clinvar_result_MEN1_2023_04_21.txt'
df = pd.read_csv(path_project + f_risultati, sep='\t', header=0)
df.drop_duplicates(inplace=True)


'''
It extracts information on whether the variation is expect to be pathogenic
 or not and dates the last review. To be checked what the date refers to.
'''

rgx_pathogenicity = r'(.+?)(?=\()'
rgx_last_review_date = r'(?<=\: )(.+?)(?=\))'

df['pathogenicity'] = df['Clinical significance (Last reviewed)'].str.extract(
    pat=rgx_pathogenicity,
    expand=True)
df['pathogenicity'] = df['pathogenicity'].astype('category')
df.pathogenicity.value_counts()

df['lr'] = df['Clinical significance (Last reviewed)'].str.extract(
    pat=rgx_last_review_date,
    expand=True)
df['lr'].dtype

df['last_review'] = pd.to_datetime(df['lr']).dt.date
df['last_review'].dtype


'''
extracts the transcript of the most relevant isoform NM_001370259.2,
which accounts for 97% of the submissions to Clinvar
'''
rgx_transcript = r'(.+?)(?=\:)'
df['transcript'] = df['Name'].str.extract(
    pat=rgx_transcript,
    expand=True)

check = df.transcript.value_counts()

rgx_parentheses = r'\s*\([^()]*\)'

df['transcript'] = df['transcript'
                      ].str.replace(rgx_parentheses,
                                    "",
                                    regex=True).str.strip()

'''
builds a dictionary of exons with their bounds
the source for the boundaries of exons is:
https://genome.ucsc.edu/cgi-bin/hgc?hgsid=1612088515_sApfRjuxpPDElNbx7W1C6iadmbkr&g=htcCdnaAliInWindow&i=NM_001370259.2&c=chr11&l=64803513&r=64810716&o=64803515&aliTable=ncbiRefSeqPsl&table=ncbiRefSeqCurated
accessed on 23 April 2023
'''
l_ranges = [
    [64810514, 64810551],
    [64810133, 64810513],
    [64809665, 64810132],
    [64808100, 64809664],
    [64807891, 64808099],
    [64807681, 64807890],
    [64807552, 64807680],
    [64807220, 64807551],    
    [64807179, 64807219],
    [64807099, 64807178],    
    [64807011, 64807098],
    [64806369, 64807010],    
    [64806232, 64806368],
    [64805771, 64806231],    
    [64805635, 64805770],
    [64805199, 64805634],    
    [64805034, 64805198],
    [64804817, 64805033],    
    [64803516, 64804816]
           ]
l_names = [
    'exon_1',
    'intron_1',
    'exon_2',
    'intron_2',
    'exon_3',
    'intron_3',
    'exon_4',
    'intron_4',
    'exon_5',
    'intron_5',
    'exon_6',
    'intron_6',
    'exon_7',
    'intron_7',
    'exon_8',
    'intron_8',
    'exon_9',
    'intron_9',
    'exon_10'
    ]

l_bins = [
    64810551,
    64810514,
    64810133,
    64809665,
    64808100,
    64807891,
    64807681,
    64807552,
    64807220,    
    64807179,
    64807099,    
    64807011,
    64806369,    
    64806232,
    64805771,    
    64805635,
    64805199,    
    64805034,
    64804817,    
    64803516
    ]
l_names.reverse()
l_bins.reverse()
#d_exons = dict(zip(l_exons, l_exon_ranges))
#d_ranges = dict(zip(l_names, l_ranges))

'''
attributes the relevant exon to each listed mutation
'''

df['var_starts'] = df.GRCh38Location.str[:9]
df['var_starts'] = df['var_starts'].astype(float).astype('Int64')
df_s = df[df.transcript == 'NM_001370259.2'].copy()

df_s.sort_values(by='var_starts', inplace=True)
df_s['exons'] = pd.cut(df_s.var_starts,
                       bins=l_bins,
                       right=False,
                       labels=l_names)


'''
extract the gene and protein variant
'''
l_pathogenic = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
check = df_s.exons.value_counts()
check_p = df_s[df_s.pathogenicity.isin(l_pathogenic)].exons.value_counts()

rgx_var = r'(?<=\:)(.+?)(?=\s*?\(|$)'

df_s['variant'] = df_s['Name'].str.extract(
    pat=rgx_var,
    expand=True)

rgx_var_protein = r'(\([^(MEN)]*\))'

df_s['protein_variant'] = df_s['Name'].str.extract(
    pat=rgx_var_protein,
    expand=True)

df_s['variant'] = df_s['Name'].str.extract(
    pat=rgx_var,
    expand=True)

'''
gets the list of single nucleotide polymorphism as per dbSNP
'''

df_s['dbSNP ID'].value_counts()
df['dbSNP ID'].unique().tolist()
s=df[df['dbSNP ID'].notnull()]['dbSNP ID'].unique().tolist()

with open('/Users/valentinaburrai/Desktop/Ids.txt', "w") as outfile:
    outfile.write("\n".join(s))

'''
loads the file with data from dbSNP, downloaded on
23 April 2023 using the list above to conduct a batch entrez extraction from
https://www.ncbi.nlm.nih.gov/sites/batchentrez
'''    

f_risults_snp = 'snp_result_NM_001370259.2_2023_04_23.txt'
df_snp = pd.read_csv(path_project + f_risults_snp, sep='\t', header=0)
df_snp.drop_duplicates(inplace=True)
df_snp = df_snp[df_snp['#chr'] != '#chr'].copy()

'''
harmonises the format of SNP ids between the two data sets
'''

df_snp['dbSNP ID'] = df_snp['snp_id'].map('rs{}'.format)

'''
merges the two datasets
'''

df_m = pd.merge(df_s, df_snp, how='outer', on='dbSNP ID')
check = df_m[df_m['dbSNP ID'].isnull()]

l_cols_to_drop = ['lr', 'Unnamed: 15', 'Clinical significance (Last reviewed)',
                  ]
df = df.drop(l_cols_to_drop, axis='columns')
