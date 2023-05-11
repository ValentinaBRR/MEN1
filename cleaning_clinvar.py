#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:55:57 2023

@author: valentinaburrai
"""
import pandas as pd
import numpy as np

'''
sets names for directories and file
'''
path_project = '/Users/valentinaburrai/Data/MEN1/Input/'
f_clinvar = 'clinvar_result_MEN1_2023_04_21.txt'
f_ensembl = 'Ensembl_exons_2023_05_01_ENST00000450708.csv'
'''
loads the file extracted from Clinvar into a pandas dataframe
'''
df = pd.read_csv(path_project + f_clinvar, sep='\t', header=0)
df.drop_duplicates(inplace=True)

'''
gets rid of all large variants affecting more than one gene,
 then discards the Gene(s) column;
gets rid of all variants that are not pathogenic

'''
df = df[df['Gene(s)'] == 'MEN1']
df = df[df['Clinical significance (Last reviewed)']
        .str.contains('pathogenic', case=False, regex=True)]
df = df[~df['Clinical significance (Last reviewed)']
        .str.contains('conflicting', case=False, regex=True)]
df.drop('Gene(s)', axis='columns', inplace=True )

'''
It extracts the transcript of the most relevant isoform NM_001370259.2,
which accounts for 97% of the submissions to Clinvar. This is considered as the
canonical transcript by Ensembl and UniProt

'''
df['name'] = df['Name']
df.name = df.name.str.replace('\(MEN1\)','', regex=True)
df['name_chromosome'] = df.name[df.name.str.contains('NC_', regex=True)]
df['name_transcript'] = df.name[df.name.str.contains('NM_', regex=True)]

rgx_transcript_name = r'(?P<transcript>N.+?(?=\:)):(?P<variant>(?<=\:).+?(?=\s*?\(|$))(?P<aa_variant>(?=\s\().+?(?<=\)))?'
df = df.join(df.name_transcript.str.extract(pat=rgx_transcript_name,
                                            expand=True))
df.aa_variant = df.aa_variant.str.strip()
df.aa_variant = df.aa_variant.str.replace('\(p.', '', regex=True
                                          ).str.replace('\)', '', regex=True)

'''
All unknown or non-canonical transcritps are dropped. 
for a list of submissions for pathogenic variants of non-canonical transcripts
see
/Users/valentinaburrai/Code/MEN1/Investigate_non_canonical_transcripts.py
'''
df = df[df.transcript == 'NM_001370259.2'].copy()

'''
separates info on pathogenicity from info on date of review
'''
rgx_pathogenicity = r'(^[^\(]+)'
rgx_last_review_date = r'(?<=\: )(.+?)(?=\))'

df['pathogenicity'] = df['Clinical significance (Last reviewed)'].str.extract(
    pat=rgx_pathogenicity,
    expand=True)
    
df['pathogenicity'] = df['pathogenicity'].str.strip()
df['pathogenicity'] = df['pathogenicity'
                         ].str.lower().str.replace(' ', '-'
                                                   ).str.replace('/', ';')

df['pathogenicity'] = df['pathogenicity'].astype('category')
df.pathogenicity.value_counts()

df['lr'] = df['Clinical significance (Last reviewed)'].str.extract(
    pat=rgx_last_review_date,
    expand=True)
df['lr'].dtype

df['last_review'] = pd.to_datetime(df['lr']).dt.date
df['last_review'].dtype

'''
drops submissions that did not specify a condition
'''

df = df[(df['Condition(s)'] != 'not provided')
        & (df['Condition(s)'] != 'not specified')]

'''
identifies the starting point of variant on gene,
Genome Reference Consortium Human Build 38
'''

df['var_starts'] = df.GRCh38Location.str[:9]
df['var_starts'] = df['var_starts'].astype(float).astype('Int64')

'''
loads the ensembl data on exon and intron bounds for the canonical transcript
http://Feb2023.archive.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000133895;r=11:64803585-64810581;t=ENST00000440873
'''
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


'''separates data referring to the main transcript'''

l_names = df_e.exon_intron[df_e.Start.notnull()].tolist()
l_bins = df_e.Start[df_e.Start.notnull()].tolist()
l_bins.append(64803516)

l_names.reverse()
l_bins.reverse()


df.sort_values(by='var_starts', inplace=True)
df['exons'] = pd.cut(df.var_starts,
                      bins=l_bins,
                      right=False,
                      labels=l_names)
           
'''
gets the list of single nucleotide polymorphism as per dbSNP


df_s['dbSNP ID'].value_counts()
df['dbSNP ID'].unique().tolist()
s=df[df['dbSNP ID'].notnull()]['dbSNP ID'].unique().tolist()

with open('/Users/valentinaburrai/Desktop/Ids.txt', "w") as outfile:
    outfile.write("\n".join(s))
'''

'''
checks that all the submissions are from the same built, so that the ch

'''
c_1 = df[~(df['Canonical SPDI'].astype(str).str.contains('NC_000011.10'))]

c_2 = df[(df['Protein change'].notnull())
           &
           (df['aa_variant'].notnull())][['Protein change', 'aa_variant']]

c_3 = df[(df['Protein change'].isnull())
             &
             (df['aa_variant'].notnull())][['Protein change', 'aa_variant']]

c_4 = df[(df['Protein change'].notnull())
             &
             (df['aa_variant'].isnull())][['Protein change', 'aa_variant']]

'''
separates the dataset in two:
    one dataset df_1 contains only one polypeptide variant,
    the other dataset df_2 lists many variants per every submission
'''
df_1 = df[df['Protein change'].notnull()].copy()
df_1 = df_1[~df_1['Protein change'].str.contains(',')].copy()

df_2 = df[df['Protein change'].notnull()].copy()
df_2 = df_2[df_2['Protein change'].str.contains(',')].copy()

df_3 = df[df['Protein change'].isna()].copy()


'''
checks that the location of the polypeptide variation coincides between name
and submission
creates one variable only with the location of the variation

'''

df_1['loc_pp_name'] = df_1.aa_variant.str.extract(r'(\d{1,3})')
df_1['loc_pp_submission'] = df_1['Protein change'].str.extract(r'(\d{1,3})')
df_1['compara'] = np.where(df_1.loc_pp_name != df_1.loc_pp_submission,
                           'problem',
                           np.nan)

df_1['aa_location'] = df_1['loc_pp_submission']

rgx_aa_1 = r'(?P<aa_original>^[A-Z]{1})(?:\d+)(?P<aa_substitution>[A-Z]{1})'
df_1 = df_1.join(df_1['Protein change'].str.extract(pat=rgx_aa_1, expand=True))


'''
explodes all variations when a submission contains more than one
creates one variable only with the location of the variation
'''
df_vars = df_2[['VariationID']].join(df_2['Protein change'].str.split(',', expand=True))
df_vars = df_vars.melt(id_vars = 'VariationID')
df_vars.drop('variable', axis='columns', inplace=True)
df_vars.rename(columns={'value': 'prot_change'}, inplace=True)
df_vars.dropna(inplace=True)
df_vars.prot_change = df_vars.prot_change.str.strip()

df_2 = pd.merge(df_2, df_vars,
                how='left',
                on='VariationID')
df_2['loc_pp_name'] = df_2.aa_variant.str.extract(r'(\d{1,3})')
df_2['loc_pp_submission'] = df_2['prot_change'].str.extract(r'(\d{1,3})')
df_2['compara'] = np.where(df_2.loc_pp_name != df_2.loc_pp_submission,
                           'problem',
                           np.nan)
check = df_2[df_2.compara != 'problem']
df_2['aa_location'] = df_2['loc_pp_submission']

#rgx_aa_2 = r'(?P<aa_original>^[A-Z]{1}[a-z]{2})(?:\d+)(?P<aa_substitution>[A-Z]{1}[a-z]{2})'
df_2 = df_2.join(df_2['prot_change'].str.extract(pat=rgx_aa_1, expand=True))

'''
extract any data on any variant
'''
df_3['aa_location'] = df_3.aa_variant.str.extract(r'(\d{1,3})')

rgx_aa_3 = r'(?P<aa_original>^[A-Z]{1}[a-z]{2})(?:\d+)(?P<aa_substitution>[A-Z]{1}[a-z]{2})'
df_3 = df_3.join(df_3['aa_variant'].str.extract(pat=rgx_aa_3, expand=True))

'''
recreates the data frame df
then etracts the position of the variant, maps it and mathces it to NextProt
'''
df= pd.concat([df_1, df_2,df_3], axis='rows', ignore_index=True)



'''
https://www.uniprot.org/uniprotkb/O00255/entry#O00255-2
maps the amino acid of the variant to the protein isoform that
NextProt considers canonical
then maps it to the functional domains listed by NextProt 
'''



#rgx_aa_var = r'(?P<aa_original>^[A-Z]{1}[a-z]{2})(?P<aa_position>\d+)(?P<aa_substitution>[A-Z]{1}[a-z]{2})'
#df = df.join(df.aa_variant.str.extract(pat=rgx_aa_var, expand=True))


df['aa_mapped_to_NextProt'] = np.where(df.aa_location.astype(float) <= 148,
                                       df.aa_location.astype(float),
                                       df.aa_location.astype(float) + 5)


df['functional_domain'] = np.where(
    df.aa_mapped_to_NextProt.between(219,395),
    'region_interacting_with_FANC2',
    np.where(df.aa_mapped_to_NextProt.between(465,502),
             'disordered_region, basic_acid_residues', np.where(
                 df.aa_mapped_to_NextProt.between(503, 557),
                 'disordered_region',
                 '')))

'''
tags each variation by type
'''
# generates a unified variable for polypeptide variants from both name and submission                            
df['aa_var'] = np.where(df.prot_change.notnull(),
                        df.prot_change,
                        df.aa_variant)   

# tags premature termination codon
df['var_type'] = np.where(
    ((df.aa_var.notnull())
     &
     ((df.aa_var.str.endswith('Ter'))
      |
      (df.aa_var.str.endswith('*')))),
                          'premature_termination_codon',
                          '')

# tags loss of termination codon
df['var_type'] = np.where(
    ((df.aa_var.notnull())
     &
     ((df.aa_var.str.startswith('Ter'))
      |
      (df.aa_var.str.startswith('*')))),
                          'loss_of_termination_codon',
                          df.var_type)

# tags synonymous variations
df['var_type'] = np.where(
    (((df.aa_var.notnull())
     &
     (df.aa_var.str.endswith('=')))
    |
    ((df.aa_substitution.notnull())
     &
     (df.aa_original == df.aa_substitution))),
                          'synonymous',
                          df.var_type)
# tags start-codon loss
df['var_type'] = np.where(
    ((df.aa_mapped_to_NextProt == 1)
     &
     (df.aa_original.isin(['Met', 'M']))),
    'start_codon_loss',
     df.var_type)

#tags splicing donor site variation
df['var_type'] = np.where(
    (
     (df.variant.str.contains(r'\+1'))
     |
     (df.variant.str.contains(r'\+2'))
     |
     (df.variant.str.contains(r'\+3'))
     ),
    'splice_donor_site_variation',
     df.var_type)

#tags splicing acceptor site variation
df['var_type'] = np.where(
    (
     (df.variant.str.contains(r'\-1'))
     |
     (df.variant.str.contains(r'\-2'))
     |
     (df.variant.str.contains(r'\-3'))
     ),
    'splice_acceptor_site_variation',
     df.var_type)

#tags frameshift variations
df['var_type'] = np.where(
    ((df.aa_var.notnull())
     &
     (df.aa_var.str.endswith('fs'))),
                          'frameshift',
                          df.var_type)

#tags remaining cases where there was a substitution of aa as missense
df['var_type'] = np.where(
    ((df.var_type == '')
     &
     (df.aa_original.notnull())
     &
     (df.aa_original != df.aa_substitution)),
                          'missense',
                          df.var_type)


d_aa_type= {
    'R': 'charged_positive',
    'H': 'charged_positive',
    'K': 'charged_positive',
    'D': 'charged_negative',
    'E': 'charged_negative',
    'S': 'polar_uncharged',
    'T': 'polar_uncharged',
    'N': 'polar_uncharged',
    'Q': 'polar_uncharged',
    'A': 'neutral',
    'V': 'neutral',
    'I': 'neutral',
    'L': 'neutral',
    'M': 'neutral',
    'F': 'neutral',
    'Y': 'neutral',
    'W': 'neutral',
    'C': 'special_case',
    'U': 'special_case',
    'G': 'special_case',
    'P': 'special_case'
    }

df['aa_original_type'] = df.aa_original.map(d_aa_type)

df['aa_substitution_type'] = df.aa_substitution.map(d_aa_type)


'''
df['aa_var'] = df.aa_variant


from functools import reduce

#str_to_replace = "The string for replacement."
#replacement_dict = {"The ": "A new ", "for ": "after "}

df['aa_var'] = reduce(lambda x, y: x.replace(*y), [df.aa_var, *list(d_aa.items())])

d_aa = {
        'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'
        }

import Bio
from Bio.Data.IUPACData import protein_letters_3to1
help(Bio.SeqUtils.IUPACData)
Bio.Data.IUPACData.protein_letters_3to1['Ala']


str_replaced = reduce(lambda x, y: x.replace(*y), [df., *list(d_aa.items())])
for k,v in d.iteritems():
    address = address.upper().replace(k, v)
'''

l_cols_to_keep = ['VariationID', 'AlleleID(s)',
                  'Review status','last_review', 'Accession',
                  'GRCh38Location', 'dbSNP ID',
                  'Canonical SPDI','pathogenicity',
                  'name', 'transcript',
                  'variant', 'var_starts', 'exons',
                  'aa_variant','prot_change',
                  'aa_original', 'aa_substitution',
                  'aa_original_type', 'aa_substitution_type',
                  'aa_location','aa_mapped_to_NextProt',
                  'functional_domain', 'var_type']

df_w = df[l_cols_to_keep].copy()

d_rename={'variant': 'gene_variant',
          'aa_variant': 'aa_variant_from_name',
          'prot_change': 'aa_variant_from_submission'}

df_w.rename(columns=d_rename, inplace=True)

df_w.sort_values('var_starts', ascending=False, inplace=True)

df_w.to_excel(path_project + 'clinvar_wrangled.xlsx', index=False,
              na_rep='',              
              freeze_panes=(1,1))
