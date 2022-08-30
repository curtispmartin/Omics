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
import json 

import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.stats.multitest as smsm # this could be a challenge... worth it? 

import matplotlib.pyplot as plt
import seaborn as sns


##### CREATE CLASS FOR HANDLING PROTEOMICS DATA
#----------------------------------------------------------------------------#
### create `dpea` class for differential protein expression analysis
class dpea:
    '''
    Handles & analyzes proteomics data. 
    
    REQUIREMENTS:
        1. Mass spectrometry (MS) data w spectral counts. 
    '''
    
### initialization parameters... SHOULD I INCLUDE THESE PARAMETERS W THE EXPERIMENT FUNCTION INSTEAD?? THINK ABOUT HOW THIS WILL BE USED (NOT AS SINGLE SCRIPT, LIKE DONE BELOW)
#     def __init__(self, df=None):
    def __init__(self, df=None, params=None):
        '''

        Parameters
        ----------
        df : TYPE, optional
            DESCRIPTION. The default is None.
        params : TYPE, optional
            DESCRIPTION. The default is params.

        Returns
        -------
        None.

        '''
        
### define some attributes which allow for data tracking
        self.source_data = df.copy()


### plot distributions for data


### process the data... remove consistently low-value entries & add pseudocounts
    def clean(self, nfloor=10, pseudocount=1):
        '''

        Parameters
        ----------
        nfloor : integer, optional
            DESCRIPTION. Threshold for removing low-value data. If the sum of all spectral counts is less than this value for this experiment (both treatments), this entry is removed. The default is 10.
        pseudocount : integer, optional
            DESCRIPTION. Value to add to spectral counts of zero in order to prevent infinite fold-change calculations. The default is 1.

        Returns
        -------
        Dataframes for each treatment.

        '''
### fill N/A w/ zeros... ASSUMES BLANK VALUES CAN BE INTERPRETED AS NOT DETECTED
        self.first_data = self.first_data.fillna(0.0)
        self.second_data = self.second_data.fillna(0.0)

### remove entries where sum of PSM count < `nfloor` (10 by default)... a heuristic which seems to work OK
#         idx_keep = ((self.first_data.sum(axis=1) >= nfloor) & (self.second_data.sum(axis=1) >= nfloor))
        idx_keep = self.first_data.join(self.second_data).sum(axis=1) >= nfloor
        self.first_data = self.first_data[idx_keep].copy()
        self.second_data = self.second_data[idx_keep].copy()

### calculate & report fraction of zeroes
        df_check = self.first_data.join(self.second_data)
        print(f'\nMissing value fraction:\n\n{(df_check == 0).sum() / df_check.shape[0]}')

### check to make sure there are no conflicts in frames
        if not self.first_data.index.equals(self.second_data.index):
            raise Exception('WARNING: Data filtering malfunctioning. Exiting...')

### add pseudocount to PSMs to prevent infinite fold-change; requires treatment of missing data, per above
        self.first_data += pseudocount
        self.second_data += pseudocount      
        
        return(self.first_data, self.second_data)
#----------------------------------------------------------------------------#
        

### normalize data (normalized spectral abundance factor, or nSAF)
    def gen_nSAF(self, plen='# AAs'):
        '''
        '''

### join protein lengths to frame 
        df1 = self.first_data.join(self.experi_data[plen]) 
        df2 = self.second_data.join(self.experi_data[plen])
    
### create new frames for normalized data
        df1_norm = pd.DataFrame()
        df2_norm = pd.DataFrame()
        
### generate normalized SAFs & keep data transformations; check to ensure values sum to one
        for col in self.first_cols:
            df1_norm[f'{col}_SAF'] = df1[col] / df1[plen] # spectral counts over protein length
            df1_norm[f'{col}_nSAF'] = df1_norm[f'{col}_SAF'] / df1_norm[f'{col}_SAF'].sum() # normalized SAF 
            if np.round(df1_norm[f'{col}_nSAF'].sum(), decimals=0) != 1:
                raise Warning(f'\nNormalized {col} does not add up to one. Check data...')

        for col in self.second_cols:            
            df2_norm[f'{col}_SAF'] = df2[col] / df2[plen] # spectral counts over protein length
            df2_norm[f'{col}_nSAF'] = df2_norm[f'{col}_SAF'] / df2_norm[f'{col}_SAF'].sum() # normalized SAF 
            if np.round(df2_norm[f'{col}_nSAF'].sum(), decimals=0) != 1:
                raise Warning(f'\nNormalized {col} does not add up to one. Check data...')
                
        self.first_data = df1[df1.columns[df1.columns != plen]].join(df1_norm)
        self.second_data = df2[df2.columns[df2.columns != plen]].join(df2_norm)
        
        return(self.first_data, self.second_data)
#----------------------------------------------------------------------------#
 

### run statistical test for determining differential expression (t-test)
    def ttest(self, alpha=0.05, correction='BH', labels=None):
        
### assign attributes for future use
        self.alpha = alpha
        self.labels = labels
                
### get nSAF columns for statistical tests & plotting later
        self.l_nsafcols1 = [col for col in self.first_data.columns if 'nSAF' in col]
        self.l_nsafcols2 = [col for col in self.second_data.columns if 'nSAF' in col]

### perform two-sample t-test using statsmodels
        ttest = stats.ttest_ind(a=self.first_data[self.l_nsafcols1], b=self.second_data[self.l_nsafcols2], alternative='two-sided', axis=1)
        l_pvals = ttest[1]

### join data for further analysis
        self.results = self.first_data.join(self.second_data)
        self.results['p-value'] = l_pvals

### implement Benjamini-Hochberg procedure to correct for multiple comparisons
#         df_sub = df_sub.sort_values('p-value') 
#         m = df_sub.shape[0] 
#         df_sub['k'] = np.arange(1, m+1, 1)
#         df_sub['p-adj'] = df_sub['p-value'] * (m / df_sub['k']) # now this should be compared against alpha
        if correction == 'BH':
            self.results['p-adj'] = smsm.multipletests(self.results['p-value'], alpha=alpha, method='fdr_bh')[1] # now this should be compared against alpha
            print(f'\nBonferroni-adjusted significance level:\t{smsm.multipletests(self.results["p-value"], alpha=alpha, method="fdr_bh")[3]:.5f}\n(just a ballpark check, Benjamini-Hochberg used here)')
        else:
            print('\nNo correction made for multiple comparisons...')

### add labels to frame if included            
        if labels:
            self.results = self.experi_data[[labels]].join(self.results, how='right')

        return(self.results)
#----------------------------------------------------------------------------#


### run statistical test for determining differential expression (LIMMA)... GONNA BE TOUGH

#----------------------------------------------------------------------------#


### generate volcano plot
    def volcano(self, fcthresh=2, alpha=None, labels=None):
        
### check for results from statistical test
        if not hasattr(self, 'results'):
            raise Exception('\nWarning: Statistical test needs to be performed before generating volcano plot. Exiting...')
            
### check to make sure certain attributes are available
#         if not hasattr(self, 'names'):
#             raise Exception('\nWarning: Missing `names` attribute. Exiting...')

### define labels if not provided explicitly in function
        if not labels:
            if hasattr(self, 'labels'):
                labels = self.labels

### new frame for plotting
        df_volc = self.experi_data[[col for col in self.experi_data.columns if col not in self.results.columns]].join(self.results, how='right')

### calculate log2(FC) = log2 of fold-change
        df_volc['Mean_Grp1'] = df_volc[self.l_nsafcols1].mean(axis=1)
        df_volc['Mean_Grp2'] = df_volc[self.l_nsafcols2].mean(axis=1)
        df_volc['FC'] = df_volc['Mean_Grp2'] / df_volc['Mean_Grp1'] # interpreting FC to mean from grp1 --> grp2

### calculate values for volcano plot (log2(FC) & log10(p-value))
        df_volc['log2(FC)'] = np.log2(df_volc['FC'])
        df_volc['-log10(p)'] = -np.log10(df_volc['p-value'])
        df_volc['-log10(p-adj)'] = -np.log10(df_volc['p-adj'])

### add flag for "statistical significance" based on user-defined alpha
        if not alpha:
            alpha = self.alpha
        df_volc.loc[df_volc[df_volc['p-value'] <= alpha].index, 'p-Sig'] = 1
        df_volc.loc[df_volc[df_volc['p-value'] > alpha].index, 'p-Sig'] = 0

## add flag for "statistical significance" based on user-defined BH-corrected alpha
        df_volc.loc[df_volc[df_volc['p-adj'] <= alpha].index, 'BH-Sig'] = 1
        df_volc.loc[df_volc[df_volc['p-adj'] > alpha].index, 'BH-Sig'] = 0

### add flag for "high fold-change" based on user-defined threshold
        fcthresh_log2 = np.log2(fcthresh)
        df_volc.loc[df_volc[np.abs(df_volc['log2(FC)']) > fcthresh_log2].index, 'HighFC'] = 1 
        df_volc.loc[df_volc[np.abs(df_volc['log2(FC)']) <= fcthresh_log2].index, 'HighFC'] = 0 

### join to rest of data, then sort for simpler selection of important genes
        df_volc = df_volc.join(self.results.loc[:, [col for col in self.results.columns if col not in df_volc.columns]])
        df_volc = df_volc.sort_values(['BH-Sig', 'p-Sig', 'HighFC', '-log10(p-adj)', 'FC'], ascending=False).reset_index(drop=True)

### add column for defining hue later
        df_volc.loc[df_volc[(df_volc['HighFC'] == 1) & (df_volc['p-Sig'] == 0)].index, 'Hue'] = 1
        df_volc.loc[df_volc[(df_volc['HighFC'] == 0) & (df_volc['p-Sig'] == 1)].index, 'Hue'] = 2
        df_volc.loc[df_volc[(df_volc['HighFC'] == 1) & (df_volc['p-Sig'] == 1)].index, 'Hue'] = 3
        df_volc.loc[df_volc[(df_volc['HighFC'] == 1) & (df_volc['BH-Sig'] == 1)].index, 'Hue'] = 4
        df_volc['Hue'] = df_volc['Hue'].fillna(0)

### set a few plot parameters... WILL NEED TO MAKE MORE FLEXIBLE
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
#         palette = {0:'#feedde', 1:'#fdbe85', 2:'#fd8d3c', 3:'#e6550d', 4:'#a63603'} # oranges
#         palette = {0:'#f2f0f7', 1:'#cbc9e2', 2:'#9e9ac8', 3:'#756bb1', 4:'#54278f'} # purples
        palette = {0:'#eff3ff', 1:'#bdd7e7', 2:'#6baed6', 3:'#3182bd', 4:'#08519c'} # blues
        markers = {0:'o', 1:'h', 2:'X', 3:'d', 4:'P'}

        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        
        sns.scatterplot(x='log2(FC)', y='-log10(p)', hue='Hue', style='Hue', data=df_volc, palette=palette, markers=markers, edgecolor='black', linewidth=0.25, s=18, legend=False, ax=ax)

        if labels:
            print('\nWorking on the labels. Can take a minute...\n')
            for idx in range(df_volc[df_volc['p-Sig'] == 1].shape[0])[:100]:
                x = 0.075+df_volc.loc[idx, 'log2(FC)']
                y = df_volc.loc[idx, '-log10(p)']
                s = df_volc.loc[idx, labels]
                plt.text(x=x, y=y, s=s, fontdict=dict(color='k', size=4))
        
        ax.set_xlim([xmin, xmax])
        ax.set_xlabel('$log_{2}$(Fold Change)')
        ax.set_ylim([ymin, ymax])
        ax.set_ylabel('$-log_{10}$(p-value)')
#         plt.title(f'Volcano Plot ({name_outp})')
        plt.title(f'Volcano Plot')
        plt.show()
        

### reorder columns & write data to file
#         if labels:
#             l_impcols = [labels] + ['FC', 'log2(FC)', 'p-value', '-log10(p)', 'p-adj', '-log10(p-adj)', 'HighFC', 'p-Sig', 'BH-Sig', 'Mean_Grp1', 'Mean_Grp2'] + self.first_cols + self.l_nsafcols1 + self.second_cols + self.l_nsafcols2
#         else:
#             l_impcols = [self.labels] + ['FC', 'log2(FC)', 'p-value', '-log10(p)', 'p-adj', '-log10(p-adj)', 'HighFC', 'p-Sig', 'BH-Sig', 'Mean_Grp1', 'Mean_Grp2'] + self.first_cols + self.l_nsafcols1 + self.second_cols + self.l_nsafcols2
#         df_volc = df_volc[l_impcols].join(df_volc.loc[:, ~df_volc.columns.isin(l_impcols)])
        
#         self.processed_data = df_volc
        
        return(fig) # return figure so user can save as s/he wants
#----------------------------------------------------------------------------#


### analyze data from experiment... CALL OTHER FUNCTIONS FOR DATA PROCESSING? ADD MORE INPUTS?
    def experiment(self, first=None, second=None, names=None): # MIGHT NEED TO INCORPORATE PLEN HERE
        '''
        Use this to do preliminary checks on data.
        - first: Columns definining first condition for testing (list-like)
        - second: Columns defining second condition for testing (list-like)
        - names: column providing protein names (str, preferably accession)
        '''
        
### assign attributes for future use
        self.names = names
        
### incorpote (& encourage!) use of json file containing experiment parameters
#         if self.params:
 
### make sure conditionsa are lists
        self.first_cols = list(first)
        self.second_cols = list(second)
        
### create new data set from source w deduplicated protein names as index
        if self.source_data[names].duplicated().sum() > 0:
            raise Warning('\nDuplicate protein names. Keeping first instance...')
            self.experi_data = self.source_data[~self.source_data.duplicated()].sort_index().copy()
        else:            
            self.experi_data = self.source_data.set_index(names).sort_index().copy()
 
### define new data split by treatment
        self.first_data = self.experi_data[self.first_cols].copy()
        self.second_data = self.experi_data[self.second_cols].copy()

### clean data & deal w missing entries
        self.clean(nfloor=10, pseudocount=1)

### normalize data
        self.gen_nSAF(plen='# AAs')

### test for statistical significance... T-TEST RIGHT NOW... ALSO DEFINITELY WANT TO BE ABLE TO USE THIS OUTSIDE OF EXPERIMENT FUNCTION
        self.ttest(alpha=0.05, correction='BH', labels='GenSymbol')
        
### generate data for overrepresentation analysis


### generate data for enrichment analysis
        
        return(self) # do I need or want this?
#----------------------------------------------------------------------------#


### pull data for testing
path_data = 'Data/KB-Appendix17-P798-16-IPs-Summary-Comparison.xlsx'
df = pd.read_excel(io=path_data, sheet_name='Summary_Comparison')

path_params = 'Data/Params/20220815_22Rv1_UN_IACS.json'
dict_params = json.load(open(path_params))

### instantiate class object using data
first = dict_params['Treatments']['1']
second = dict_params['Treatments']['2']
names = 'Accession'
test = dpea(df=df)
test.experiment(first=first, second=second, names=names)

df_first = test.first_data
df_second = test.second_data
df_results = test.results

fig = test.volcano()
fig.savefig('test.png', bbox_inches='tight', dpi=300)


sys.exit()


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
    
    

#----------------------------------------------------------------------------#


##### LOAD, CLEAN, & SELECT DATA
#----------------------------------------------------------------------------#


### define & make output directory... helps keep things organized!
path_outp = os.path.join(os.getcwd(), 'Output', name_inpu)
if not os.path.exists(path_outp):
    os.mkdir(path_outp)


#----------------------------------------------------------------------------#



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



