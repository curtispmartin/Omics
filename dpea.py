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
import statsmodels.stats.multitest as smsm # only to check adjusted p-values, for now

import matplotlib.pyplot as plt
import seaborn as sns


##### MAXIMIZE USER INPUT TO IMPROVE WORKFLOW GENERALIZATION
#----------------------------------------------------------------------------#
### get & format user inputs
filepath = os.path.abspath(sys.argv[1])
dict_args = dict(arg.split('=') for arg in sys.argv[2:])
# print(dict_args)

### set confidence for determining "statistical significance"
if '--alpha' in dict_args.keys():
    alpha = float(dict_args['--alpha'])
    print(f'\nAlpha = {alpha:0.2f}...')
else:
    alpha = 0.05
    print(f'\nSetting alpha to {alpha:0.2f}...')


### NEW FEATURE: IMPORT PARAMETER FILE FOR SMOOTHER PROCESSING


### set threshold for determining "high fold-change"
if '--fc' in dict_args.keys():
    fcthresh = float(dict_args['--fc'])
    print(f'Fold-change = {fcthresh:0.2f}...')
else:
    fcthresh = 2.0
    print(f'Setting fold-change threshold to {fcthresh:0.2f}...')

### set flag for labeling points (0 means no label; 1 means label)
if '--label' in dict_args.keys():
    label_flag = 0
    print('Not labeling points...')
else:
    label_flag = 1

### get filename & path to data from user input
# datafiles = pd.Series(os.listdir(path_data))
# filename = datafiles.loc[int(input(f'\n{datafiles}\n\nSelect index of input file for analysis: '))].split('.')[0]
filename = os.path.split(filepath)[-1].split('.')[0]
path_data = os.path.split(filepath)[0]

### have the user select which sheet in data
# datasheets = pd.Series(pd.ExcelFile(f'{path_data}/{filename}.xlsx').sheet_names)
datasheets = pd.Series(pd.ExcelFile(filepath).sheet_names)
sheetname = datasheets.loc[int(input(f'\n{datasheets}\n\nSelect index of sheet to use for analysis:\t'))]
#----------------------------------------------------------------------------#


##### LOAD, CLEAN, & SELECT DATA
#----------------------------------------------------------------------------#
### load data
# data = pd.read_excel(io=f'{path_data}/{filename}.xlsx', sheet_name=sheetname)
data = pd.read_excel(io=filepath, sheet_name=sheetname)

### define & make output directory... helps keep things organized!
# path_outp = f'{os.getcwd()}/Output/{filename}'
path_outp = os.path.join(os.getcwd(), 'Output', filename)
if not os.path.exists(path_outp):
    os.mkdir(path_outp)

### identify condition columns either through script call or via interactive input
if '--cond1' in dict_args.keys():
    l_psmcols1 = dict_args['--cond1'].split(',')
else:
    l_psmcols1 = input(f'\n{pd.Series(data.columns)}\n\nSelect PSM data for first sample:\t\t').split()

if '--cond2' in dict_args.keys():
    l_psmcols2 = dict_args['--cond2'].split(',')
else:
    l_psmcols2 = input(f'\n{pd.Series(data.columns)}\n\nSelect PSM data for second sample:\t\t').split()

### get desired name for output files
oname = input('Provide name for output files:\t\t\t')

### convert indices to integers
l_psmcols1 = [int(i) for i in l_psmcols1]
l_psmcols2 = [int(i) for i in l_psmcols2]

### convert indices to columns
l_psmcols1 = data.columns[l_psmcols1]
l_psmcols2 = data.columns[l_psmcols2]

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
# print(df1.describe())
# print(df2.describe())
#----------------------------------------------------------------------------#


##### RUN STATISTICAL TEST, CORRECTED FOR MULTIPLE HYPOTHESIS TESTING
#----------------------------------------------------------------------------#
### create function for generating normalize spectral abundance factor (nSAF)
def normalize(df, data, col, plen='# AAs'):

### join amino acid counts to frame for metric generation
#     df = df.join(data['# AAs']) # THIS IS NOT NECESSARILY A STANDARD COLUMN!! NEED TO GENERALIZE!
    df = df.join(data[plen]) # THIS IS NOT NECESSARILY A STANDARD COLUMN!! NEED TO GENERALIZE!

### generate new metrics needed to normalize data
    df_norm = pd.DataFrame()
    df_norm[f'{col}_SAF'] = df[col] / df['# AAs'] # spectral counts over protein length
    df_norm[f'{col}_nSAF'] = df_norm[f'{col}_SAF'] / df_norm[f'{col}_SAF'].sum() # normalized SAF 
#     df_norm[f'{col}_nSAF'] = df_norm[f'{col}_SAF'] / df_norm[f'{col}_SAF'].max() # THIS IS WRONG!!!
#     print(df_norm[f'{col}_nSAF'].sum())
    
    return(df_norm)

### normalize PSM --> nSAF
for col in df1.columns:
    df1 = df1.join(normalize(df1, data, col))
for col in df2.columns:
    df2 = df2.join(normalize(df2, data, col))
    
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
# df_sub = df_sub.sort_index()

print(f'\nBonferroni-adjusted significance level:\t{smsm.multipletests(df_sub["p-value"], alpha=alpha, method="fdr_bh")[3]:.5f}\n(just a ballpark check, Benjamini-Hochberg used here)')
#----------------------------------------------------------------------------#


##### PREPARE DATA FOR PLOTTING & OUTPUT
#----------------------------------------------------------------------------#
### filter data to mitochondrial proteins only... MAY NEED TO MAKE MORE FLEXIBLE, DEPENDS ON DATA CONSISTENCY & RELEVANCE
# f_mito = 'Mitochondrial' 
# df_vol = data[data['Mitochondrial location'] == f_mito].copy()
df_vol = data.copy()

### create new frame for volcano plot
df_vol = df_vol[['Accession','GenSymbol']].join(df_sub)

### calculate log2(FC) = log2 of fold-change
df_vol['Mean_Grp1'] = df_vol[l_nsafcols1].mean(axis=1)
df_vol['Mean_Grp2'] = df_vol[l_nsafcols2].mean(axis=1)
df_vol['FC'] = df_vol['Mean_Grp2'] / df_vol['Mean_Grp1'] # interpreting FC to mean from grp1 --> grp2

### calculate values for volcano plot (log2(FC) & log10(p-value))
df_vol['log2(FC)'] = np.log2(df_vol['FC'])
df_vol['-log10(p)'] = -np.log10(df_vol['p-value'])
df_vol['-log10(p-adj)'] = -np.log10(df_vol['p-adj'])

### add flag for "statistical significance" based on user-defined alpha
# df_vol.loc[df_vol[df_vol['p-value'] <= alpha].index, 'Significance'] = 'Y'
# df_vol.loc[df_vol[df_vol['p-value'] > alpha].index, 'Significance'] = 'N'

### add flag for "high fold-change" based on user-defined threshold
fcthresh_log2 = np.log2(fcthresh)
df_vol.loc[df_vol[np.abs(df_vol['log2(FC)']) > fcthresh_log2].index, 'HighFC'] = 1 
df_vol.loc[df_vol[np.abs(df_vol['log2(FC)']) <= fcthresh_log2].index, 'HighFC'] = 0 

### join to rest of data, then sort for simpler selection of important genes
df_vol = df_vol.join(data.loc[:, [col for col in data.columns if col not in df_vol.columns]])
# df_vol = df_vol.sort_values(['HighFC', 'FC', '-log10(p-adj)'], ascending=False).reset_index(drop=True)
df_vol = df_vol.sort_values(['HighFC', '-log10(p-adj)', 'FC'], ascending=False).reset_index(drop=True)
#----------------------------------------------------------------------------#


##### PLOT & SAVE DATA TO FILE
#----------------------------------------------------------------------------#
### set a few plot parameters
xmin = -5
xmax = 5
ymin = 0
ymax = np.max(df_vol['-log10(p)']) + 0.1

### generate plot parameters based on user-defined parameters above
xmin_alpha = 0.20
xmax_alpha = fcthresh_log2 / (xmax - xmin)
ymin_fc = -np.log10(alpha) / (ymax - ymin)
ymax_fc = 0.95

### make a volcano plot w/ p-values
palette = {1:'firebrick', 0:'Silver'}
fig, ax = plt.subplots(1, 1, figsize=(8,6))

sns.scatterplot(x='log2(FC)', y='-log10(p)', hue='HighFC', data=df_vol, palette=palette, s=20, edgecolor=None, lw=0, ax=ax)
# ax.axhline(y=-np.log10(alpha), xmin=xmin_alpha, xmax=0.5-xmax_alpha, color='C3', ls='--')
# ax.axhline(y=-np.log10(alpha), xmin=0.5+xmax_alpha, xmax=1.0-xmin_alpha, color='C3', ls='--')
# ax.axvline(x=-fcthresh_log2, ymin=ymin_fc, ymax=ymax_fc, color='C3', ls='--')
# ax.axvline(x=fcthresh_log2, ymin=ymin_fc, ymax=ymax_fc, color='C3', ls='--')

if label_flag == 1:
    print('\nWorking on the labels...\n')
    for idx in range(df_vol.shape[0])[:100]:
        x = 0.05+df_vol.loc[idx, 'log2(FC)']
        y = df_vol.loc[idx, '-log10(p)']
        s = df_vol.loc[idx, 'GenSymbol']
        if (pd.isnull(x)) | (pd.isnull(y)) | (pd.isnull(s)):
            pass
        else:
            plt.text(x=x, y=y, s=s, fontdict=dict(color='k', size=5))

ax.set_xlim([xmin, xmax])
ax.set_xlabel('$log_{2}$(Fold Change)')
ax.set_ylim([ymin, ymax])
ax.set_ylabel('$-log_{10}$(p-value)')
# plt.legend(title=f'$log_{2}$(FC) > {fcthresh_log2:0.2f}', loc='upper right')
plt.legend(title=f'Fold change > {fcthresh}', loc='lower right', ncol=2)
plt.title(f'Volcano Plot ({oname})')

if label_flag == 1:
    plt.savefig(f'{path_outp}/{oname}_volcano.png', bbox_inches='tight', dpi=300)
else:
    plt.savefig(f'{path_outp}/{oname}_volcano-nolabel.png', bbox_inches='tight', dpi=300)
# plt.show()

### reorder columns & write data to file
if label_flag == 1:
    l_impcols = ['Accession', 'GenSymbol', 'FC', 'log2(FC)', 'p-value', '-log10(p)', 'p-adj', '-log10(p-adj)', 'Protein names']
    df_vol = df_vol[l_impcols].join(df_vol.loc[:, ~df_vol.columns.isin(l_impcols)])
    df_vol.to_excel(f'{path_outp}/{oname}_processed.xlsx', index=False)


### make a volcano plot w/ adjusted p-values
ymax = np.max(df_vol['-log10(p-adj)']) + 0.1

# palette = {1:'firebrick', 0:'Silver'}
fig, ax = plt.subplots(1, 1, figsize=(8,6))

sns.scatterplot(x='log2(FC)', y='-log10(p-adj)', hue='HighFC', data=df_vol, palette=palette, s=20, edgecolor=None, lw=0, ax=ax)
# ax.axhline(y=-np.log10(alpha), xmin=xmin_alpha, xmax=0.5-xmax_alpha, color='C3', ls='--')
# ax.axhline(y=-np.log10(alpha), xmin=0.5+xmax_alpha, xmax=1.0-xmin_alpha, color='C3', ls='--')
# ax.axvline(x=-fcthresh_log2, ymin=ymin_fc, ymax=ymax_fc, color='C3', ls='--')
# ax.axvline(x=fcthresh_log2, ymin=ymin_fc, ymax=ymax_fc, color='C3', ls='--')

if label_flag == 1:
#     print('\nWorking on the labels...\n')
    for idx in range(df_vol.shape[0])[:100]:
        x = 0.05+df_vol.loc[idx, 'log2(FC)']
        y = df_vol.loc[idx, '-log10(p-adj)']
        s = df_vol.loc[idx, 'GenSymbol']
        if (pd.isnull(x)) | (pd.isnull(y)) | (pd.isnull(s)):
            pass
        else:
            plt.text(x=x, y=y, s=s, fontdict=dict(color='k', size=5))

ax.set_xlim([xmin, xmax])
ax.set_xlabel('$log_{2}$(Fold Change)')
ax.set_ylim([ymin, ymax])
ax.set_ylabel('$-log_{10}$(Adjusted p-value)')
# plt.legend(title=f'$log_{2}$(FC) > {fcthresh_log2:0.2f}', loc='upper right')
plt.legend(title=f'Fold change > {fcthresh}', loc='lower right', ncol=2)
plt.title(f'Volcano Plot ({oname})')

if label_flag == 1:
    plt.savefig(f'{path_outp}/{oname}_volcano-adj.png', bbox_inches='tight', dpi=300)
else:
    plt.savefig(f'{path_outp}/{oname}_volcano-adj-nolabel.png', bbox_inches='tight', dpi=300)
# plt.show()

### create output frame
df_outp = df_vol.sort_values('FC').copy()

### provide reference list for performing analyses against
df_outp['Accession'].to_csv(f'{path_outp}/{oname}_ref.txt', sep='\t', header=False, index=False)

### provide ranked list of proteins for use w overrepresentation analysis (ORA) in PANTHER
df_outp[df_outp['HighFC'] == 1]['Accession'].to_csv(f'{path_outp}/{oname}_ora.txt', sep='\t', header=False, index=False)

### provide ranked list of proteins for use w enrichment (FCA) analysis in PANTHER
df_outp[['Accession', 'FC']].dropna().to_csv(f'{path_outp}/{oname}_sea.txt', sep='\t', header=False, index=False)
# df_outp.to_string(f'{path_outp}/{oname}_ranked.txt', header=False, index=False)
# with open(f'{path_outp}/{oname}_ranked.txt', 'w') as text:
#     for idx, entry in df_outp.iterrows():
#         entryA = entry['Accession']
#         entryB = entry['FC']
#         text.write(f'{entryA}\t{entryB}\n')
#----------------------------------------------------------------------------#  



