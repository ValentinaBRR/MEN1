#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:55:57 2023

@author: valentinaburrai
"""
import pandas as pd

'''
sets names for directories and file
'''
path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
f_clinvar = 'clinvar_result_MEN1_2023_04_21.txt'

'''
loads the file extracted from Clinvar into a pandas dataframe
'''
df = pd.read_csv(path_project + f_clinvar, sep='\t', header=0)
df.drop_duplicates(inplace=True)

'''
It extracts information on whether the variation is expect to be pathogenic
 or not and dates the last review from the clinvar variable 
 Clinical significance (Last reviewed).
To be checked what the date refers to.
'''

rgx_pathogenicity = r'(.+?)(?=\()'
rgx_last_review_date = r'(?<=\: )(.+?)(?=\))'

df['pathogenicity'] = df['Clinical significance (Last reviewed)'].str.extract(
    pat=rgx_pathogenicity,
    expand=True)

d_pathogenicity = {
    'Uncertain significance': 'uncertain',
     'Likely benign': 'likely_benign',
     'Pathogenic': 'pathogenic',
     'Conflicting interpretation of pathogenicity': 'conflicting_pathogenicity',
     'Likely pathogenic': 'likely_pathogenic',
     'Benign': 'benign',
     'Pathogenic/Likely pathogenic': 'pathogenic_likely_pathogenic',
     'Benign/Likely benign': 'benign_likely_benign',
                   }
    
df['pathogenicity'] = df['pathogenicity'].map(d_pathogenicity)
df['pathogenicity'] = df['pathogenicity'].astype('category')
df.pathogenicity.value_counts()

df['lr'] = df['Clinical significance (Last reviewed)'].str.extract(
    pat=rgx_last_review_date,
    expand=True)
df['lr'].dtype

df['last_review'] = pd.to_datetime(df['lr']).dt.date
df['last_review'].dtype


'''
It extracts the transcript of the most relevant isoform NM_001370259.2,
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
identifies the starting point of variant
'''

df['var_starts'] = df.GRCh38Location.str[:9]
df['var_starts'] = df['var_starts'].astype(float).astype('Int64')

'''
extract the gene and protein variant
'''

rgx_var = r'(?<=\:)(.+?)(?=\s*?\(|$)'

df['variant'] = df['Name'].str.extract(
   pat=rgx_var,
   expand=True)

rgx_var_protein = r'(\([^(MEN)]*\))'

df['protein_variant'] = df['Name'].str.extract(
   pat=rgx_var_protein,
   expand=True)

df['variant'] = df['Name'].str.extract(
   pat=rgx_var,
   expand=True)

check = df[(df.protein_variant.notnull())
           &
           (df['Protein change'].notnull())]

check = df[(df.protein_variant.notnull())
           &
           (df['Protein change'].isnull())]

'''separates data referring to the main transcript'''

'''
builds lsits of exons and intron with their bounds
the source for the boundaries of exons is:
https://genome.ucsc.edu/cgi-bin/hgc?hgsid=1612088515_sApfRjuxpPDElNbx7W1C6iadmbkr&g=htcCdnaAliInWindow&i=NM_001370259.2&c=chr11&l=64803513&r=64810716&o=64803515&aliTable=ncbiRefSeqPsl&table=ncbiRefSeqCurated
accessed on 23 April 2023
'''
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


df_s = df[df.transcript == 'NM_001370259.2'].copy()

df_s.sort_values(by='var_starts', inplace=True)
df_s['exons'] = pd.cut(df_s.var_starts,
                       bins=l_bins,
                       right=False,
                       labels=l_names)

l_pathogenic = ['pathogenic', 'pathogenic_likely_pathogenic',
                'likely_pathogenic']

check = df_s.exons.value_counts()
check_p = df_s[df_s.pathogenicity.isin(l_pathogenic)].exons.value_counts()
           


'''
gets the list of single nucleotide polymorphism as per dbSNP


df_s['dbSNP ID'].value_counts()
df['dbSNP ID'].unique().tolist()
s=df[df['dbSNP ID'].notnull()]['dbSNP ID'].unique().tolist()

with open('/Users/valentinaburrai/Desktop/Ids.txt', "w") as outfile:
    outfile.write("\n".join(s))
'''

l_cols_to_drop = ['lr', 'Unnamed: 15', 'Clinical significance (Last reviewed)']
df_s.drop(l_cols_to_drop, axis='columns', inplace=True)

df_s.to_excel(path_project + 'clinvar_wrangled.xlsx', index=False)
