#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 10:46:40 2023

@author: valentinaburrai
"""
import pandas as pd

'''
sets name for directory
'''
path_project = '/Users/valentinaburrai/Data/MEN1/Input/'


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

d_ranges = dict(zip(l_names, l_ranges))
df = pd.DataFrame.from_dict(d_ranges, orient='index', columns=['lower_bound',
                                                               'higher_bound'])
df.to_excel(path_project + 'NM_001370259.2Exons.xlsx')
