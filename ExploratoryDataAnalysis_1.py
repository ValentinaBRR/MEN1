#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 13 19:20:23 2023

@author: valentinaburrai
"""

import pandas as pd
import numpy as np

'''
sets names for directories and file
'''
path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
f_clinvar = 'cv_full.xlsx'

'''
loads Clinvar data tagged with Ensembl, UniProt and NextProt information
and tagged for type of variant
'''
df = pd.read_excel(path_project + f_clinvar)
df.drop_duplicates(inplace=True)


df_summary = df.describe(include='all', datetime_is_numeric=True)


df_multiple_submissions = df[df.duplicated(subset='var_id', keep=False)]

df_polypeptide = df[(df.aa_var.notnull())
                    &
                    (~ df.var_type.isin(['splice_acceptor_site_variation', 'splice_donor_site_variation']))].copy()

df_same_aa_var_from_multiple_submissions = df[df.duplicated(
    subset=['aa_location', 'aa_var'], keep=False)]

df_exons_by_domain = df[['functional_domain','exons']                    
                         ].groupby(['functional_domain',
                                    'exons']
                                    )['exons'].count()

df_pathogenic = df.groupby('pathogenicity')['var_id'].count()
check = df[df.pathogenicity == 'pathogenic;likely-pathogenic']

df_evidence = df.groupby('evidence_status')['var_id'].count()
df_evidence = pd.pivot_table()

check = df[df.pathogenicity == 'pathogenic;likely-pathogenic']

df_confidence_1 = pd.crosstab(index=df.exons,
                            columns=df.pathogenicity,
                            margins=True,
                            normalize=True)
                                    

df_no_polypeptide = df[df.aa_var.isnull()]

df_no_poly_no_intron = df[df.aa_var.isnull()
                          &
                          df.var_type.isnull()]            

df['evidence_status'].value_counts()

df_same_aa_var_from_multiple_submissions['evidence_status'].value_counts()

df.missense_type.value_counts()
