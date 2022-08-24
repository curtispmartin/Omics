#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 12:22:56 2022

@author: Curtis P. Martin

Parameters
--------------
alpha: p-value considered "statistically significant". Currently only used for plotting
fcthresh: fold-change considered "signficant". log2() transformation applied later

"""

### import modules
import os
import sys

import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as smsm # to calculate BH-adjusted p-values
import json # for reading experiment parameters

import matplotlib.pyplot as plt
import seaborn as sns


##### MAXIMIZE USER INPUT TO IMPROVE WORKFLOW GENERALIZATION
#----------------------------------------------------------------------------#
dict_args = dict(arg.split('=') for arg in sys.argv[1:])


### NEW FEATURE: IMPORT PARAMETER FILE FOR SMOOTHER PROCESSING
if '-i' in dict_args.keys():
    dict_params = json.load(open(dict_args['-i']))

### set confidence level for "statistical significance"
    alpha = float(dict_params['Parameters']['alpha'])
#     print(f'\nAlpha = {alpha:0.2f}...')
    
### set threshold for determining "high fold-change"
    fcthresh = float(dict_params['Parameters']['fcthresh'])
#     print(f'Fold-change = {fcthresh:0.2f}...')
    
### set flag for labeling points (0 means no label; 1 means label)
    label_flag = int(dict_params['Parameters']['label'])

### get columns required for analysis
    l_psmcols1 = dict_params['Treatments']['1']
    l_psmcols2 = dict_params['Treatments']['2']

### get column for normalizing PSMs against protein length
    plen = dict_params['Metadata']['Lengths']
    
### get IO variables
    name_inpu = os.path.split(dict_params['Metadata']['Data'])[-1].split('.')[0]
    name_sheet = dict_params['Metadata']['Sheet']
    path_data = os.path.split(dict_params['Metadata']['Data'])[0]
    path_inpu = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0], dict_params['Metadata']['Data'])
    name_outp = dict_params['Metadata']['Name']

### load data
    data = pd.read_excel(io=path_inpu, sheet_name=name_sheet)
    
    
### otherwise provide inputs via bash
else:
    path_inpu = os.path.abspath(dict_args['-f'])
    
### set confidence level for "statistical significance"
    if '--alpha' in dict_args.keys():
        alpha = float(dict_args['--alpha'])
    else:
        alpha = 0.05
    
### set threshold for determining "high fold-change"
    if '--fc' in dict_args.keys():
        fcthresh = float(dict_args['--fc'])
    else:
        fcthresh = 2.0
    
### set flag for labeling points (0 means no label; 1 means label)
    if '--label' in dict_args.keys():
        label_flag = 0
        print('Not labeling points...')
    else:
        label_flag = 1

### get filename & path to data from user input
    name_inpu = os.path.split(path_inpu)[-1].split('.')[0]
    path_data = os.path.split(path_inpu)[0]

### have the user select which sheet in data
    datasheets = pd.Series(pd.ExcelFile(path_inpu).sheet_names)
    name_sheet = datasheets.loc[int(input(f'\n{datasheets}\n\nSelect index of sheet to use for analysis:\t'))]
    
### identify condition columns either through script call or via interactive input
    if '--cond1' in dict_args.keys():
        l_psmcols1 = dict_args['--cond1'].split(',')
    else:
        l_psmcols1 = input(f'\n{pd.Series(data.columns)}\n\nSelect PSM data for first sample:\t\t').split()
    
    if '--cond2' in dict_args.keys():
        l_psmcols2 = dict_args['--cond2'].split(',')
    else:
        l_psmcols2 = input(f'\n{pd.Series(data.columns)}\n\nSelect PSM data for second sample:\t\t').split()

### load data
    data = pd.read_excel(io=path_inpu, sheet_name=name_sheet)
        
### convert indices to integers
    l_psmcols1 = [int(i) for i in l_psmcols1]
    l_psmcols2 = [int(i) for i in l_psmcols2]
    
### convert indices to columns
    l_psmcols1 = data.columns[l_psmcols1]
    l_psmcols2 = data.columns[l_psmcols2]
    
### get desired name for output files
    name_outp = input('Provide name for output files:\t\t\t')
#----------------------------------------------------------------------------#


##### LOAD, CLEAN, & SELECT DATA
#----------------------------------------------------------------------------#
### provide feedback to user
print(f'\nSetting alpha to:\t\t\t{alpha:0.2f}')
print(f'Setting fold-change threshold to:\t{fcthresh:0.2f}')

### define & make output directory... helps keep things organized!
path_outp = os.path.join(os.getcwd(), 'Output', name_inpu)
if not os.path.exists(path_outp):
    os.mkdir(path_outp)

### select data of interest
df1 = data[l_psmcols1].copy()
df2 = data[l_psmcols2].copy()

### fill N/A w/ zeros... ASSUMES BLANK VALUES CAN BE INTERPRETED AS NOT DETECTED
df1 = df1.fillna(0.0)
df2 = df2.fillna(0.0)

### remove entries where sum of PSM count < 10... SEEMS TO BE A DECENT COMPROMISE HEURISTIC
nfloor = 10
# idx_keep = ((df1.sum(axis=1) >= nfloor) & (df2.sum(axis=1) >= nfloor))
idx_keep = df1.join(df2).sum(axis=1) >= nfloor
df1 = df1[idx_keep].copy()
df2 = df2[idx_keep].copy()

### remove entries where all values are zero... IT THIS REASONABLE?
# print(df1.join(df2).loc[(~(df1==0).all(axis=1)) & (~(df2==0).all(axis=1))])
# idx_keep = (~(df1==0).all(axis=1)) & (~(df2==0).all(axis=1))
# df1 = df1[idx_keep].copy()
# df2 = df2[idx_keep].copy()

### calculate & report fraction of zeroes
df_check = df1.join(df2)
print(f'\nFraction of zero-values found:\n\n{(df_check == 0).sum() / df_check.shape[0]}')

### check to make sure there are no conflicts in frames
if not df1.index.equals(df2.index):
    raise Exception('WARNING: Data filtering malfunctioning. Exiting...')

### add pseudo-count to PSMs to remove infinite fold-change; requires treatment of missing data, per above
df1 += 1
df2 += 1
#----------------------------------------------------------------------------#


##### RUN STATISTICAL TEST, CORRECTED FOR MULTIPLE HYPOTHESIS TESTING
#----------------------------------------------------------------------------#
### create function for generating normalize spectral abundance factor (nSAF)
def normalize(df, data, col, plen='# AAs'):

### join amino acid counts to frame for metric generation
    df = df.join(data[plen]) # WILL NEED TO GENERALIZE NORMALIZATION CONSTANT INPUT

### generate new metrics needed to normalize data
    df_norm = pd.DataFrame()
    df_norm[f'{col}_SAF'] = df[col] / df['# AAs'] # spectral counts over protein length
    df_norm[f'{col}_nSAF'] = df_norm[f'{col}_SAF'] / df_norm[f'{col}_SAF'].sum() # normalized SAF 
    
    return(df_norm)

### normalize PSM --> nSAF
for col in df1.columns:
    df1 = df1.join(normalize(df1, data, col, plen=plen))
for col in df2.columns:
    df2 = df2.join(normalize(df2, data, col, plen=plen))
    
### get nSAF columns for statistical tests & plotting later
l_nsafcols1 = [col for col in df1.columns if 'nSAF' in col]
l_nsafcols2 = [col for col in df2.columns if 'nSAF' in col]

### perform two-sample t-test on PSM data
ttest = stats.ttest_ind(a=df1[l_nsafcols1], b=df2[l_nsafcols2], alternative='two-sided', axis=1)
# print(ttest[1])
# l_tstats = ttest[0]
l_pvals = ttest[1]

### perform Mann-Whitney U test on PSM data
mwutest = stats.mannwhitneyu(x=df1[l_nsafcols1], y=df2[l_nsafcols2], alternative='two-sided', axis=1)
# print(mwutest[1])
# l_mwustats = mwutest[0]
# l_pvals = mwutest[1]

### join data for further analysis
df_sub = df1.join(df2)
df_sub['p-value'] = l_pvals

### implement Benjamini-Hochberg procedure to correct for multiple comparisons
# df_sub = df_sub.sort_values('p-value') 
# m = df_sub.shape[0] 
# df_sub['k'] = np.arange(1, m+1, 1)
# df_sub['p-adj'] = df_sub['p-value'] * (m / df_sub['k']) # now this should be compared against alpha
df_sub['p-adj'] = smsm.multipletests(df_sub['p-value'], alpha=alpha, method='fdr_bh')[1] # now this should be compared against alpha
print(f'\nBonferroni-adjusted significance level:\t{smsm.multipletests(df_sub["p-value"], alpha=alpha, method="fdr_bh")[3]:.5f}\n(just a ballpark check, Benjamini-Hochberg used here)')
#----------------------------------------------------------------------------#


##### PREPARE DATA FOR PLOTTING & OUTPUT
#----------------------------------------------------------------------------#
### filter data to mitochondrial proteins only... MAY NEED TO MAKE MORE FLEXIBLE, DEPENDS ON DATA CONSISTENCY & RELEVANCE
# f_mito = 'Mitochondrial' 
# df_volc = data[data['Mitochondrial location'] == f_mito].copy()
df_volc = data.copy()

### create new frame for volcano plot
df_volc = df_volc[['Accession','GenSymbol']].join(df_sub)

### calculate log2(FC) = log2 of fold-change
df_volc['Mean_Grp1'] = df_volc[l_nsafcols1].mean(axis=1)
df_volc['Mean_Grp2'] = df_volc[l_nsafcols2].mean(axis=1)
df_volc['FC'] = df_volc['Mean_Grp2'] / df_volc['Mean_Grp1'] # interpreting FC to mean from grp1 --> grp2

### calculate values for volcano plot (log2(FC) & log10(p-value))
df_volc['log2(FC)'] = np.log2(df_volc['FC'])
df_volc['-log10(p)'] = -np.log10(df_volc['p-value'])
df_volc['-log10(p-adj)'] = -np.log10(df_volc['p-adj'])

### add flag for "statistical significance" based on user-defined alpha
df_volc.loc[df_volc[df_volc['p-value'] <= alpha].index, 'p-Sig'] = 1
df_volc.loc[df_volc[df_volc['p-value'] > alpha].index, 'p-Sig'] = 0

## add flag for "statistical significance" based on user-defined BH-corrected alpha
df_volc.loc[df_volc[df_volc['p-adj'] <= alpha].index, 'BH-Sig'] = 1
df_volc.loc[df_volc[df_volc['p-adj'] > alpha].index, 'BH-Sig'] = 0

### add flag for "high fold-change" based on user-defined threshold
fcthresh_log2 = np.log2(fcthresh)
df_volc.loc[df_volc[np.abs(df_volc['log2(FC)']) > fcthresh_log2].index, 'HighFC'] = 1 
df_volc.loc[df_volc[np.abs(df_volc['log2(FC)']) <= fcthresh_log2].index, 'HighFC'] = 0 

### add flag for high fold-change & high significance, based on user definitions
# df_volc.loc[df_volc[df_volc[['HighFC', 'p-Sig']].sum(axis=1) == 2].index, 'POI'] = 1
# df_volc.loc[df_volc[df_volc[['HighFC', 'p-Sig']].sum(axis=1) != 2].index, 'POI'] = 0

### join to rest of data, then sort for simpler selection of important genes
df_volc = df_volc.join(data.loc[:, [col for col in data.columns if col not in df_volc.columns]])
df_volc = df_volc.sort_values(['BH-Sig', 'p-Sig', 'HighFC', '-log10(p-adj)', 'FC'], ascending=False).reset_index(drop=True)

### add column for defining hue later
df_volc.loc[df_volc[(df_volc['HighFC'] == 1) & (df_volc['p-Sig'] == 0)].index, 'Hue'] = 1
df_volc.loc[df_volc[(df_volc['HighFC'] == 0) & (df_volc['p-Sig'] == 1)].index, 'Hue'] = 2
df_volc.loc[df_volc[(df_volc['HighFC'] == 1) & (df_volc['p-Sig'] == 1)].index, 'Hue'] = 3
df_volc.loc[df_volc[(df_volc['HighFC'] == 1) & (df_volc['BH-Sig'] == 1)].index, 'Hue'] = 4
df_volc['Hue'] = df_volc['Hue'].fillna(0)
#----------------------------------------------------------------------------#


##### PLOT & SAVE DATA TO FILE
#----------------------------------------------------------------------------#
### set a few plot parameters
xmin = -5
xmax = 5
ymin = 0
ymax = np.max(df_volc['-log10(p)']) + 0.1

### generate plot parameters based on user-defined parameters above
xmin_alpha = 0.20
xmax_alpha = fcthresh_log2 / (xmax - xmin)
ymin_fc = -np.log10(alpha) / (ymax - ymin)
ymax_fc = 0.95

### make a volcano plot w/ p-values
# palette = {0:'#feedde', 1:'#fdbe85', 2:'#fd8d3c', 3:'#e6550d', 4:'#a63603'} # oranges
# palette = {0:'#f2f0f7', 1:'#cbc9e2', 2:'#9e9ac8', 3:'#756bb1', 4:'#54278f'} # purples
palette = {0:'#eff3ff', 1:'#bdd7e7', 2:'#6baed6', 3:'#3182bd', 4:'#08519c'} # blues
markers = {0:'o', 1:'h', 2:'X', 3:'d', 4:'P'}
fig, ax = plt.subplots(1, 1, figsize=(8,6))

sns.scatterplot(x='log2(FC)', y='-log10(p)', hue='Hue', style='Hue', data=df_volc, palette=palette, markers=markers, edgecolor='black', linewidth=0.25, s=18, legend=False, ax=ax)
if label_flag == 1:
    print('\nWorking on the labels. Can take a minute...\n')
    for idx in range(df_volc[df_volc['p-Sig'] == 1].shape[0])[:100]:
        x = 0.075+df_volc.loc[idx, 'log2(FC)']
        y = df_volc.loc[idx, '-log10(p)']
        s = df_volc.loc[idx, 'GenSymbol']
        plt.text(x=x, y=y, s=s, fontdict=dict(color='k', size=4))

ax.set_xlim([xmin, xmax])
ax.set_xlabel('$log_{2}$(Fold Change)')
ax.set_ylim([ymin, ymax])
ax.set_ylabel('$-log_{10}$(p-value)')
plt.title(f'Volcano Plot ({name_outp})')

if label_flag == 1:
    plt.savefig(f'{path_outp}/{name_outp}_volcano.png', bbox_inches='tight', dpi=300)
else:
    plt.savefig(f'{path_outp}/{name_outp}_volcano-nolabel.png', bbox_inches='tight', dpi=300)
# plt.show()

### reorder columns & write data to file
if label_flag == 1:
    l_impcols = ['Accession', 'GenSymbol', 'FC', 'log2(FC)', 'p-value', '-log10(p)', 'p-adj', '-log10(p-adj)', 'Protein names', 'HighFC', 'p-Sig', 'BH-Sig']
    df_volc = df_volc[l_impcols].join(df_volc.loc[:, ~df_volc.columns.isin(l_impcols)])
    df_volc.to_excel(f'{path_outp}/{name_outp}_processed.xlsx', index=False)

### make a volcano plot w/ adjusted p-values
# ymax = np.max(df_volc['-log10(p-adj)']) + 0.1

# palette = {1:'firebrick', 0:'Silver'}
# fig, ax = plt.subplots(1, 1, figsize=(8,6))

# sns.scatterplot(x='log2(FC)', y='-log10(p-adj)', hue='Hue', style='Hue', data=df_volc, palette=palette, markers=markers, edgecolor='black', linewidth=0.25, s=15, legend=False, ax=ax)
# if label_flag == 1:
#     print('\nWorking on the labels. They can take a minute...\n')
#     for idx in range(df_volc[df_volc['p-Sig'] == 1].shape[0])[:100]:
#         x = 0.075+df_volc.loc[idx, 'log2(FC)']
#         y = df_volc.loc[idx, '-log10(p)']
#         s = df_volc.loc[idx, 'GenSymbol']
#         plt.text(x=x, y=y, s=s, fontdict=dict(color='k', size=4))

# ax.set_xlim([xmin, xmax])
# ax.set_xlabel('$log_{2}$(Fold Change)')
# ax.set_ylim([ymin, ymax])
# ax.set_ylabel('$-log_{10}$(p-value)')
# plt.title(f'Volcano Plot ({name_outp})')

# ax.set_xlim([xmin, xmax])
# ax.set_xlabel('$log_{2}$(Fold Change)')
# ax.set_ylim([ymin, ymax])
# ax.set_ylabel('$-log_{10}$(Adjusted p-value)')
# plt.title(f'Volcano Plot ({name_outp})')

# if label_flag == 1:
#     plt.savefig(f'{path_outp}/{name_outp}_volcano-adj.png', bbox_inches='tight', dpi=300)
# else:
#     plt.savefig(f'{path_outp}/{name_outp}_volcano-adj-nolabel.png', bbox_inches='tight', dpi=300)
# plt.show()

### create frame for printing to output
df_outp = df_volc.sort_values('FC').copy()

### get list of all proteins detected in experiment for use as reference in further analyses
# df_outp['Accession'].to_csv(f'{path_outp}/{name_outp}_ref.txt', sep='\t', header=False, index=False)
df_outp_ref = data[(data[l_psmcols1 + l_psmcols2].fillna(0.0).sum(axis=1) != 0)]['Accession'].sort_values() # can't use df1/df2 b/c of filtering procedure... need all proteins detected in experiment
# print(df_outp_ref)
df_outp_ref.to_csv(f'{path_outp}/{name_outp}_ref.txt', sep='\t', header=False, index=False)

### provide ranked list of proteins for use w overrepresentation analysis (ORA) in PANTHER
df_outp[df_outp['HighFC'] == 1]['Accession'].to_csv(f'{path_outp}/{name_outp}_ora.txt', sep='\t', header=False, index=False)

### provide ranked list of proteins for use w enrichment (FCA) analysis in PANTHER
df_outp[['Accession', 'FC']].dropna().to_csv(f'{path_outp}/{name_outp}_sea.txt', sep='\t', header=False, index=False)
#----------------------------------------------------------------------------#  



