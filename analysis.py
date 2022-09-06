#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Thu Aug 18 12:27:25 2022

@author: martincup

'''

### import modules
import os
import sys
import json 

import pandas as pd
import numpy as np

import scipy.stats as stats
# import scipy.interpolate as scint # for fitting natural cubic spline in calculating q-values
import sklearn.isotonic as skliso # utilize the isotonic regression function in place of natural cubic spline which isn't available
import statsmodels.stats.multitest as smsm # this could be a challenge... worth it? 

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/Users/martincup/Research/Omics') # add Omics directory to python path
import omics # work in progress!

import requests # for accessing uniprot api


### define experimental parameters
path_params = 'Params/20220815_22Rv1_UN_IACS.json' # the more interesting data 
# path_params = 'Params/20220815_PC3_UN_006.json'
dict_params = json.load(open(path_params))

### get IO variables
name_inpu = os.path.split(dict_params['Metadata']['Data'])[-1].split('.')[0]
name_sheet = dict_params['Metadata']['Sheet']
path_data = os.path.split(dict_params['Metadata']['Data'])[0]
path_inpu = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0], dict_params['Metadata']['Data'])
name_outp = dict_params['Metadata']['Name']
path_outp = os.path.join(os.getcwd(), 'Output', name_inpu, name_outp)

### create output directory if none exists
if not os.path.exists(path_outp):
    os.makedirs(path_outp)


### get data used for enrichment analysis
# names = ['accession', 'fold-change']
# df_sea = pd.read_csv(f'Output/KB-Appendix17-P798-16-IPs-Summary-Comparison/sea-{name_outp}.txt', sep='\t', names=names, header=None)

### get processed data generated from `dpea` pull
df_proc = pd.read_csv(os.path.join(path_outp, f'processed-{name_outp}.csv'))

### get enriched pathway data
names = ['PANTHER Pathways', '#', '+/-', 'P value', 'FDR']
df_enri = pd.read_csv(os.path.join(path_outp, f'enriched-{name_outp}.txt'), sep='\t', names=names, header=4)

### get pathway accession codes
l_enri = df_enri['PANTHER Pathways'].str.split().str[-1].str.split('(').str[-1].str.split(')').str[0].tolist()
df_enri['Pathway Accession'] = l_enri
df_enri = df_enri[['Pathway Accession'] + names].copy()


### get gene list derived from enrichment data via PANTHER
names = ['Pathway Accession', 'Mapped IDs', 'Pathway Name', 'Components', 'Subfamilies', 'Associated']
df_list = pd.read_csv(os.path.join(path_outp, f'genelist-{name_outp}.txt'), sep='\t', names=names, header=None)

### break down list into separate protein identifiers (using 'Mapped IDs')
dict_prot = {}
for idx, row in df_list.iterrows():
    str_ids = row['Mapped IDs'].split(',')
    l_ids = []
    for item in str_ids:
        l_ids.append(item.split('|')[-1].split('=')[-1])    
    dict_prot[row['Pathway Accession']] = l_ids
    
    
### get proteins associated w enriched pathways
df_prot = pd.DataFrame()
for idx, row in df_enri.iterrows():
    pathway = row['Pathway Accession']
    description = row['PANTHER Pathways']

    df_int = df_proc[df_proc['Accession'].isin(dict_prot[pathway])].copy()
    df_int['Pathway'] = pathway
    df_int['Description'] = description
    df_int = df_int.sort_values(by='GenSymbol')
    
    df_prot = pd.concat([df_prot, df_int])


### define function for pulling protein name & information from UniProt API
def get_protname(accession='P60709', organism_id=9606):
    
### define url for pulling protein name from API... IS THERE AN EASIER WAY? 
    url = f'https://rest.uniprot.org/uniprotkb/search?query=accession:{accession}&organism_id:{organism_id}&columns=protein_name'    

### request uniprot API & print status
    response = requests.get(url)
    if not response.status_code == 200:
        raise Exception(f'\nCaution: Link to UniProt API failed for {accession}. Exiting...')
    
### pull data of interest... right now, just uniprotid (for confirmation) & protein name
    uniprotid = response.json()['results'][0]['uniProtkbId']
    protname = response.json()['results'][0]['proteinDescription']['recommendedName']['fullName']['value']

### message for user
    print(f'Working on {accession}... {protname}')
    
### close request... not sure this is strictly necessary
    response.close()
    
    return(protname)

### test it out
# accession = 'P60709'
# organism_id = 9606 # human
# print(get_protname(accession=accession, organism_id=organism_id))

### get protein names for proteins in enriched sets
df_prot['Protein Name'] = df_prot['Accession'].apply(lambda accession: get_protname(accession=accession, organism_id=9606))

### reorganize for simple viewing
df_prot = df_prot[['Pathway', 'Description', 'Accession', 'GenSymbol', 'fold-change', 'p-value', 'q-value', 'Protein Name']].reset_index(drop=True).copy()

### save to file
df_prot.to_csv(os.path.join(path_outp, f'protlist-{name_outp}.csv'), index=False)


# sys.exit()


### format data for dot plot
# fcthresh = 2.0
df_plot = df_prot.sort_values('fold-change', ascending=False)
df_plot['size'] = 1 / df_plot['q-value']
df_plot.loc[df_plot[df_plot['fold-change'] >= 1.0].index, 'direction'] = '+'
df_plot.loc[df_plot[df_plot['fold-change'] <= 1.0].index, 'direction'] = '-'
df_plot.loc[df_plot[df_plot['direction'] == '-'].index, 'fold-change'] = 1 / df_plot['fold-change'] # essentially want to plot absolute value for easier interpretation
nplot = 50
df_plot = pd.concat([df_plot[df_plot['direction'] == '+'].head(n=nplot), df_plot[df_plot['direction'] == '-'].tail(n=nplot)])
df_plot = df_plot.sort_values(by='GenSymbol')
df_plot['Description'] = df_plot['Description'].str.split().str[:2].str.join(sep=' ')
# df_plot = df_plot.sort_values(by='fold-change', ascending=False)

### a hybrid b/w a heat map & a dot plot
sns.set_theme(style='whitegrid')
g = sns.relplot(data=df_plot, x='GenSymbol', y='Description', hue='q-value', size='fold-change', palette='Reds_r', hue_norm=(0, 1), edgecolor='k', height=10, sizes=(50, 250), size_norm=(1, 10))

### improve plot formatting
g.set(xlabel='', ylabel='', title='Proteins in Enriched Pathways', aspect='equal')
g.despine(left=True, bottom=True)
g.ax.margins(.05)
for label in g.ax.get_xticklabels():
    label.set_rotation(90)
for artist in g.legend.legendHandles:
    artist.set_edgecolor('k')

g.savefig(os.path.join(path_outp, f'heatmap-{name_outp}.png'), dpi=300)


sys.exit()


dict_count = df_plot['GenSymbol'].value_counts().to_dict()
for idx, row in df_plot.iterrows():
    df_plot.loc[idx, 'Count'] = dict_count[row['GenSymbol']]

g = sns.relplot(x='q-value', y='fold-change', hue='Pathway', size='Count', data=df_plot)
g.savefig(os.path.join(path_outp, f'test-{name_outp}.png'), dpi=300)


sys.exit()


### create a radial dot plot? 
sns.set_theme()
g = sns.FacetGrid(df_plot, hue='Pathway', subplot_kws=dict(projection='polar'), height=4.5, sharex=False, sharey=False, despine=False)
g.map(sns.scatterplot, 'q-value', 'fold-change')
# g.savefig(os.path.join(path_outp, f'test-{name_outp}.png'), dpi=300)










 