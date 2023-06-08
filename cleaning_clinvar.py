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

rgx_transcript_name = r'(?P<transcript>N.+?(?=\:)):(?P<gene_var>(?<=\:).+?(?=\s*?\(|$))(?P<aa_var>(?=\s\().+?(?<=\)))?'
df = df.join(df.name_transcript.str.extract(pat=rgx_transcript_name,
                                            expand=True))
df.aa_var = df.aa_var.str.strip()
df.aa_var = df.aa_var.str.replace('\(p.', '', regex=True
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
df['last_review'] = pd.to_datetime(df['lr'])


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

df['aa_location'] = df.aa_var.str.extract(r'(\d{1,3})')

rgx_aa_var = r'(?P<aa_wild_type>^[A-Z]{1}[a-z]{2})(?:\d+)(?P<aa_mutated>[A-Z]{1}[a-z]{2})'
df = df.join(df['aa_var'].str.extract(pat=rgx_aa_var, expand=True))


'''
https://www.uniprot.org/uniprotkb/O00255/entry#O00255-2
maps the amino acid of the variant to the protein isoform that
NextProt considers canonical
then maps it to the functional domains listed by NextProt 
'''


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
#tags delins
df['var_type'] = np.where(
    ((df.aa_var.notnull())
     &
     (df.aa_var.str.contains('del'))),
                          'delins',
                          '')

#tags small deletes
df['var_type'] = np.where(
    ((df.aa_var.notnull())
     &
     (df.aa_var.str.contains('del'))),
                          'deletion',
                          df.var_type)

# tags premature termination codon
df['var_type'] = np.where(
    ((df.aa_var.notnull())
     &
     ((df.aa_var.str.endswith('Ter'))
      |
      (df.aa_var.str.endswith('*')))),
                          'premature_termination_codon',
                          df.var_type)

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
    ((df.aa_mutated.notnull())
     &
     (df.aa_wild_type == df.aa_mutated))),
                          'synonymous',
                          df.var_type)
# tags start-codon loss
df['var_type'] = np.where(
    ((df.aa_mapped_to_NextProt == 1)
     &
     (df.aa_wild_type.isin(['Met', 'M']))),
    'start_codon_loss',
     df.var_type)

#tags splicing donor site variation
df['var_type'] = np.where(
     (df.gene_var.str.contains(r'\+\d{1}')),
    'splice_donor_site_variant',
     df.var_type)

#tags splicing acceptor site variation
df['var_type'] = np.where(
     (df.gene_var.str.contains(r'\-\d{1}')),
    'splice_acceptor_site_variant',
     df.var_type)

#tags intronic variant
df['var_type'] = np.where(
     ((df.gene_var.str.contains(r'\-\d{2}'))
      |(df.gene_var.str.contains(r'\+\d{2}'))),
    'intronic_variant',
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
     (df.aa_wild_type.notnull())
     &
     (df.aa_wild_type != df.aa_mutated)),
                          'missense',
                          df.var_type)

'''
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
    'F': 'aromatic',
    'Y': 'aromatic',
    'W': 'aromatic',
    'C': 'polar_uncharged',
    'G': 'neutral',
    'P': 'neutral'
    }
'''
#{k:v for k, v in d_aa.items() if v == 'P'}
d_aa_type= {
    'Arg': 'charged_positive',
    'His': 'charged_positive',
    'Lys': 'charged_positive',
    'Asp': 'charged_negative',
    'Glu': 'charged_negative',
    'Ser': 'polar_uncharged',
    'Thr': 'polar_uncharged',
    'Asn': 'polar_uncharged',
    'Gln': 'polar_uncharged',
    'Ala': 'neutral',
    'Val': 'neutral',
    'Ile': 'neutral',
    'Leu': 'neutral',
    'Met': 'neutral',
    'Phe': 'aromatic',
    'Tyr': 'aromatic',
    'Trp': 'aromatic',
    'Cys': 'polar_uncharged',
    'Gly': 'neutral',
    'Pro': 'neutral'
    }

df['aa_wt_type'] = df.aa_wild_type.map(d_aa_type)

df['aa_m_type'] = df.aa_mutated.map(d_aa_type)

df['missense_type'] = np.where(
    (df.aa_wt_type.notnull()
    &
    df.aa_m_type.notnull()),
    df.aa_wt_type + '_to_' + df.aa_m_type,
    ''
    )


d_rename={'VariationID': 'var_id',
          'Condition(s)': 'condition',
          'AlleleID(s)': 'allele_id',
          'dbSNP ID': 'dbSNP_id',
          'Canonical SPDI': 'canonical_SPDI',
          'Accession': 'accession',
          'GRCh38Location': 'GRCh38_location',
          }

df.rename(columns=d_rename, inplace=True)

'''
d_aa = {
        'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'
        }
'''
'''
finds duplicates among the polypeptide variations which were submitted in
different submissions.
It then ranks them by the quality of the evidence and when the quality of the
evidence is missing it ranks them by submission id.
It keeps the submission with the highest evidence rank, when ranks are equal
it keeps those submitted first.
'''
d_rev_status = {
    'no assertion provided': 
        'The allele was included in a submission that did not provide'
        ' an interpretation.',
    'no assertion criteria provided': 
        'The allele was included in a submission with an interpretation'
        ' but without assertion criteria.',
    'no assertion for the individual variant':
        'The allele was not interpreted  directly in any submission; '
        'it was submitted to ClinVar only as a component of a '
        'compound heterozygote or a haplotype.',
    'criteria provided, single submitter': 
        'One star: One submitter provided an interpretation with assertion'
        ' criteria.',
    'criteria provided, conflicting interpretations':
        'One star: Multiple Submitters provided assertion criteria but there'
        ' are conflicting interpretations. The independent values are'
        ' enumerated for clininal significance.',
    'criteria provided, multiple submitters, no conflicts':
        'Two stars: Two or more submitters with assertion criteria provided'
        ' the same interpretation.',
    'reviewed by expert panel':
        'Three stars: the variant was reviewed by an expert panel.',
    'Practice guideline':
        'Four stars: The variant was reviewed by a professional society that'
        ' provides practice guideliines.'
                 }
    
evidence_status = pd.Categorical(df['Review status'], categories=d_rev_status.keys(),
                                 ordered=True)
df['evidence_status'] = evidence_status 

df_check = df[(df['Review status'] != df.evidence_status)]

df_same_aa_var_from_multiple_submissions = df[df.duplicated(
    subset=['aa_location', 'aa_var'], keep=False)].copy()

df_same_aa_var_from_multiple_submissions.sort_values(
    by=['aa_location', 'aa_var', 'evidence_status', 'accession'],
    ascending=[True, True, False, True],
    inplace=True)

df_keep = df_same_aa_var_from_multiple_submissions.drop_duplicates(
    subset=['aa_location', 'aa_var'], keep='first')
l_to_drop = df_same_aa_var_from_multiple_submissions.index[
    ~df_same_aa_var_from_multiple_submissions.index.isin(df_keep.index)]


'''
generates three outputs: 
    one output for genes (as in some instances there are more than one gene
                          variant per protein variant)
    one output for proteins (dropping duplicate genes variants for the same
                             protein variant)
    one output for introni variants
                             
'''
l_cols_to_keep = [
    'var_id', 'allele_id',
    'name','dbSNP_id',
    'canonical_SPDI',
    'last_review', 'accession',
    'GRCh38_location', 'pathogenicity',
    'condition',
    'evidence_status',
    'gene_var', 'var_starts', 'exons',
    'aa_location','aa_mapped_to_NextProt',
    'aa_var',
    'aa_wild_type', 'aa_mutated',
    'aa_wt_type', 'aa_m_type',
    'var_type',
    'functional_domain',
    'missense_type'
    ]
df_full = df[l_cols_to_keep].copy()
df_full.sort_values('var_starts', ascending=False, inplace=True)
df_full.to_excel(path_project + 'cv_full.xlsx', index=False,
                                 na_rep='',              
                                 freeze_panes=(1,1))



l_cols_to_keep_intronic = [
    'var_id', 'allele_id',
    'name','dbSNP_id',
    'canonical_SPDI',
    'last_review', 'accession',
    'GRCh38_location', 'pathogenicity',
    'evidence_status',
    'gene_var', 'var_starts', 'exons',
    'var_type', 'missense_type'
                            ]

df_intr = df[df.var_type.str.contains('splic')][l_cols_to_keep_intronic]
df_intr.sort_values('var_starts', ascending=False, inplace=True)

df_intr.to_excel(path_project + 'cv_intronic.xlsx', index=False,
              na_rep='',              
              freeze_panes=(1,1))

l_cols_to_keep_gene = [
    'var_id', 'allele_id',
    'name','dbSNP_id',
    'canonical_SPDI',
    'last_review', 'accession',
    'GRCh38_location', 'pathogenicity',
    'evidence_status',
    'gene_var', 'var_starts', 'exons',
                        ]
df_gene = df[~df.var_type.str.contains('splic')][l_cols_to_keep_gene]
df_gene.sort_values('var_starts', ascending=False, inplace=True)

df_gene.to_excel(path_project + 'cv_gene.xlsx', index=False,
              na_rep='',              
              freeze_panes=(1,1))

l_cols_to_keep_protein = [
    'var_id', 'allele_id',
    'name','dbSNP_id',
    'canonical_SPDI',
    'last_review', 'accession',
    'GRCh38_location', 'pathogenicity',
    'evidence_status', 'exons',
    'aa_var',
    'aa_wild_type', 'aa_mutated',
    'aa_wt_type', 'aa_m_type',
    'aa_location','aa_mapped_to_NextProt',
    'functional_domain', 'var_type', 'missense_type'
                            ]

df_prot = df[~df.index.isin(l_to_drop)][l_cols_to_keep_protein]
df_prot.sort_values('aa_mapped_to_NextProt', ascending=True, inplace=True)
df_prot.to_excel(path_project + 'cv_protein.xlsx', index=False,
              na_rep='',              
              freeze_panes=(1,1))
