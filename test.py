#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Aug 29 22:43:22 2022

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

# sys.path.append('/Users/martincup/Research/Omics') # add Omics directory to python path
path_work = os.getcwd()
sys.path.append(path_work) # add Omics directory to python path
import omics # work in progress!

import requests # for accessing APIs


### pull data for testing
path_data = 'Data/KB-Appendix17-P798-16-IPs-Summary-Comparison.xlsx'
df = pd.read_excel(io=path_data, sheet_name='Summary_Comparison')

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


### open file for writing run commentary
# path_comments = os.path.join(path_outp, 'output.txt')


##### DIFFERENTIAL EXPRESSION ANALYSIS
#----------------------------------------------------------------------------#
### create dpea object
exp = omics.dpea(df=df)

### define experiment
first = dict_params['Treatments']['1']
second = dict_params['Treatments']['2']
names = 'Accession'
exp = exp.experiment(first=first, second=second, names=names)

### clean data
exp = exp.clean(nfloor=10, pseudocount=1)

### normalize data
plen = dict_params['Metadata']['Lengths']
exp = exp.norm_nSAF(plen=plen)

### run statistical test
# log = True
alpha = float(dict_params['Parameters']['alpha']) # confidence level
exp = exp.ttest(alpha=alpha, correction='BH', labels='GenSymbol')

### plot test data distributions
fig = exp.plot_distr()
fig.savefig(os.path.join(path_outp, f'distr-{name_outp}.png'), bbox_inches='tight', dpi=300)

### try running all steps in compound method
# exp = omics.dpea(df=df).experiment(first=first, second=second, names=names).clean(nfloor=10, pseudocount=1).norm_nSAF(plen='# AAs').ttest(alpha=0.05, correction='BH', labels='GenSymbol')

### generate volcano plot
fcthresh = float(dict_params['Parameters']['fcthresh']) # threshold for relevant fold-change
fig = exp.volcano(fcthresh=fcthresh, labels='GenSymbol')
fig.savefig(os.path.join(path_outp, f'volcano-{name_outp}.png'), bbox_inches='tight', dpi=300)

### generate q-values & some diagnostic plots
exp = exp.est_qval(pval='p-value')

fig = exp.plot_pvq()
# fig.savefig(os.path.join(path_outp, f'pvq-{name_outp}.png'), bbox_inches='tight', dpi=300)

fig = exp.plot_sigvq()
# fig.savefig(os.path.join(path_outp, f'sigvq-{name_outp}.png'), bbox_inches='tight', dpi=300)

### prepare data for ORA
df_ora, df_ref = exp.prep_ora()
df_ora.to_csv(os.path.join(path_outp, f'ora-{name_outp}.txt'), index=False, header=False)
df_ref.to_csv(os.path.join(path_outp, f'ref-{name_outp}.txt'), index=False, header=False)

### prepare data for SEA
df_sea = exp.prep_sea()
df_sea.to_csv(os.path.join(path_outp, f'sea-{name_outp}.txt'), sep='\t', index=False, header=False)

### save data for fidelity
df_outp = exp.results.copy()
df_outp.loc[df_outp[df_outp['q-value'] < alpha].index, 'significance'] = 1
df_outp['significance'] = df_outp['significance'].fillna(0)
df_outp = df_outp.sort_values(by=['significance', 'fold-change'], ascending=False)
df_outp.to_csv(os.path.join(path_outp, f'processed-{name_outp}.csv'))
#----------------------------------------------------------------------------#


##### STATISTICAL ENRICHMENT ANALYSIS (SEA) USING PANTHERDB API
#----------------------------------------------------------------------------#
### get some parameters for analysis
organism = 9606 # human
correction = 'FDR' # p-value correction via false discovery rate
annot = 'ANNOT_TYPE_ID_PANTHER_PATHWAY' # pathway set to search (e.g., Panther or GO Biological Process)... CODE THESE BETTER
path_data = os.path.join(path_outp, f'sea-{name_outp}.txt')
    
### run enrichment analysis & save to file
enri = omics.enrichment(path_data=path_data)
df_enri = enri.run_sea(cutoff=0.05)
df_enri.to_csv(os.path.join(path_outp, f'enriched-{name_outp}.csv'), index=False)

### match proteins in data to pathways identified by enrichement analysis
df_protid = enri.match_sea()
df_protid.to_csv(os.path.join(path_outp, f'genelist-{name_outp}.csv'), index=False)

### get processed data generated from `dpea` pull... CAN I GET RID OF THIS NOW THAT IT'S IN THE SAME SCRIPT?
df_proc = pd.read_csv(os.path.join(path_outp, f'processed-{name_outp}.csv'))

### get enriched pathway data... CAN I GET RID OF THIS NOW THAT IT'S IN THE SAME SCRIPT?
df_enri = pd.read_csv(os.path.join(path_outp, f'enriched-{name_outp}.csv'))

### get gene list derived from enrichment data via PANTHER... CAN I GET RID OF THIS NOW THAT IT'S IN THE SAME SCRIPT?
df_list = pd.read_csv(os.path.join(path_outp, f'genelist-{name_outp}.csv'))

### merge data
# df_prot = df_proc.merge(df_list.merge(df_enri, on='Pathway'), on='Accession', suffixes=(None, '_mergecheck'), indicator=True).sort_values(by=['FDR', 'q-value'])
df_prot = df_proc.merge(df_list.merge(df_enri, on='Pathway'), on='Accession', suffixes=(None, '_mergecheck')).sort_values(by=['FDR', 'q-value'])

### get protein names for proteins in enriched sets
df_prot['Protein Name'] = df_prot['Accession'].apply(lambda accession: enri.get_protname(accession=accession, organism_id=9606))

### save results
df_prot.to_csv(os.path.join(path_outp, f'results-{name_outp}.csv'), index=False)
#----------------------------------------------------------------------------#


sys.exit()


##### DOT PLOT FOR VISUALIZATION OF ENRICHMENT RESULTS
#----------------------------------------------------------------------------#
### format data for dot plot
df_plot = df_prot.sort_values('fold-change', ascending=False)
df_plot.loc[df_plot[df_plot['fold-change'] >= 1.0].index, 'direction'] = '+'
df_plot.loc[df_plot[df_plot['fold-change'] <= 1.0].index, 'direction'] = '-'
df_plot.loc[df_plot[df_plot['direction'] == '-'].index, 'fold-change'] = -1 / df_plot['fold-change'] # essentially want to plot absolute value for easier interpretation
nplot = 50
df_plot = pd.concat([df_plot[df_plot['direction'] == '+'].head(n=nplot), df_plot[df_plot['direction'] == '-'].tail(n=nplot)])
df_plot = df_plot.sort_values(by='GenSymbol')
df_plot['Description'] = df_plot['Description'].str.split().str[:2].str.join(sep=' ')
df_plot['-log10q'] = -np.log10(df_plot['q-value'])

### a hybrid b/w a heat map & a dot plot
sns.set_theme(style='whitegrid')
# g = sns.relplot(data=df_plot, x='GenSymbol', y='Description', hue='q-value', size='fold-change', palette='Reds_r', hue_norm=(0, 1), edgecolor='k', height=10, sizes=(50, 250), size_norm=(1, 10))
g = sns.relplot(data=df_plot, x='GenSymbol', y='Description', hue='fold-change', size='-log10q', palette='vlag', hue_norm=(-5, 5), edgecolor='k', height=10, sizes=(50, 250), size_norm=(0, 2))

### improve plot formatting
g.set(xlabel='', ylabel='', title='Proteins in Enriched Pathways', aspect='equal')
g.despine(left=True, bottom=True)
g.ax.margins(.05)
for label in g.ax.get_xticklabels():
    label.set_rotation(90)
for artist in g.legend.legendHandles:
    artist.set_edgecolor('k')

g.savefig(os.path.join(path_outp, f'dot-{name_outp}.png'), dpi=300)
#----------------------------------------------------------------------------#






sys.exit()


### implement linear model for microarray data (LIMMA) method for analyzing differential expression

### create dpea object
test = omics.dpea(df=df)

### define experiment
first = dict_params['Treatments']['1']
second = dict_params['Treatments']['2']
names = 'Accession'
test = test.experiment(first=first, second=second, names=names)


sys.exit()


### format data for dot plot
df_plot = df_outp[['GenSymbol', 'fold-change', 'q-value']].sort_values('fold-change', ascending=False)
df_plot['size'] = 1 / df_plot['q-value']
df_plot.loc[df_plot[df_plot['q-value'] < alpha].index, 'color'] = 'C1'
df_plot['color'] = df_plot['color'].fillna('C0')
df_plot.loc[df_plot[df_plot['fold-change'] >= fcthresh].index, 'direction'] = '+'
df_plot.loc[df_plot[df_plot['fold-change'] <= fcthresh].index, 'direction'] = '-'
df_plot.loc[df_plot[df_plot['direction'] == '-'].index, 'fold-change'] = 1 / df_plot['fold-change'] # essentially want to plot absolute value for easier interpretation

nplot = 50
df_plot = pd.concat([df_plot[df_plot['direction'] == '+'].head(n=nplot), df_plot[df_plot['direction'] == '-'].tail(n=nplot)])

### create dot plot using relplot function
g = sns.relplot(x='fold-change', y='GenSymbol', hue='color', palette={'C0':'C0', 'C1':'C1'}, size='size', col='direction', data=df_plot, edgecolor='k', height=10, aspect=.5, legend=False, facet_kws={'sharey': False, 'sharex': True})

### use same x axis limits on all columns and add better labels
g.set(xlim=[0,np.round(np.max(df_plot['fold-change']) + 1)], xlabel='Fold-Change', ylabel='')

### use more meaningful titles for the columns
titles = ['Upregulated', 'Downregulated']

for ax, title in zip(g.axes.flat, titles):

### set a different title for each axis
    ax.set(title=title)

### make the grid horizontal instead of vertical
#     ax.xaxis.grid(False)
#     ax.yaxis.grid(True, lw=1, color='lightgrey')

sns.despine(left=True, bottom=True)

g.savefig(os.path.join(path_outp, f'dot-{name_outp}.png'), bbox_inches='tight', dpi=300)










    

