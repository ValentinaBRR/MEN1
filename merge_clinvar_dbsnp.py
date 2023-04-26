#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:55:57 2023

@author: valentinaburrai
"""

import pandas as pd


'''
sets names for directory and files
'''
path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
f_risults_snp = 'snp_result_NM_001370259.2_2023_04_23.txt'
f_c_wrangled = 'clinvar_wrangled.xlsx'

'''
loads wrangled clinvar
'''
df_c = pd.read_excel(path_project + f_c_wrangled)

'''
loads the file with data from dbSNP, downloaded on
23 April 2023 using the list above to conduct a batch entrez extraction from
https://www.ncbi.nlm.nih.gov/sites/batchentrez
see cleaning_clinvar.py
'''    

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

df_m = pd.merge(df_c, df_snp, how='outer', on='dbSNP ID')
check = df_m[df_m['dbSNP ID'].isnull()]

df_m.to_excel(path_project+'clinvar_dbsnp.xlsx', index=False)
