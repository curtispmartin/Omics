#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Tue Apr 12 12:22:56 2022

@author: Curtis P. Martin

'''

### import modules
import os
import sys
import json 

import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.stats.multitest as smsm # this could be a challenge... worth it? 
import sklearn.isotonic as skliso # utilize the isotonic regression function in place of natural cubic spline which isn't available

import matplotlib.pyplot as plt
import seaborn as sns

import requests # for accessing APIs


##### CREATE CLASS FOR HANDLING PROTEOMICS DATA
#----------------------------------------------------------------------------#
### create `dpea` class for differential protein expression analysis
class dpea:
    '''
    Handles & analyzes proteomics data. 
    
    REQUIREMENTS:
        1. Mass spectrometry (MS) data w spectral counts. 
    '''
    
#----------------------------------------------------------------------------#
### initialization parameters... SHOULD I INCLUDE THESE PARAMETERS W THE EXPERIMENT FUNCTION INSTEAD?? THINK ABOUT HOW THIS WILL BE USED (NOT AS SINGLE SCRIPT, LIKE DONE BELOW)
    def __init__(self, df=None):
#     def __init__(self, df=None, params=None):
        '''

        Parameters
        ----------
        df : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        '''
        
### define some attributes which allow for data tracking
        self.source_data = df.copy()
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### set everything up for interpreting & analyzing experimental data
    def experiment(self, first=None, second=None, names=None): 
        '''
        
        Parameters
        ----------
        first : list, required
            List of columns definining first condition for testing. The default is None.
        second : list, required
            List of columns defining second condition for testing. The default is None.
        names : str, required
            Name for column defining some sort of protein identifier (e.g., accession). The default is None.

        Returns
        -------
        self.

        '''

### assign attributes for future use
        self.names = names
 
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
        
        return(self) # return self to enable compound method calls
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### process the data... remove consistently low-value entries & add pseudocounts
    def clean(self, nfloor=10, pseudocount=1, **kwargs):     
        '''

        Parameters
        ----------
        nfloor : integer, optional
            Threshold for removing low-value data. If the sum of all spectral counts is less than this value for this experiment (both treatments), this entry is removed. The default is 10.
        pseudocount : integer, optional
            Value to add to spectral counts of zero in order to prevent infinite fold-change calculations. The default is 1.
        **kwargs : pd.DataFrame, optional
            Can provide new data for cleaning using 'first' & 'second' options, like in experiment() method. 

        Returns
        -------
        self.

        '''
        
### check to see if user provided new data for cleaning
        if 'first' in kwargs.keys():
            if hasattr(self, 'first_data'):
                print('\nWarning: Copying over data provided in experiment() call.')
            self.first_data = kwargs['first']            
            
        if 'second' in kwargs.keys():
            if hasattr(self, 'second_data'):
                print('\nWarning: Copying over data provided in experiment() call.')
            self.second_data = kwargs['second']

### fill N/A w/ zeros... ASSUMES BLANK VALUES CAN BE INTERPRETED AS NOT DETECTED
        self.first_data = self.first_data.fillna(0.0)
        self.second_data = self.second_data.fillna(0.0)

### remove entries where sum of PSM count < `nfloor` (10 by default)... a heuristic which seems to work OK
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
        
        return(self)
#----------------------------------------------------------------------------#
        

#----------------------------------------------------------------------------#
### normalize data (normalized spectral abundance factor, or nSAF)
    def norm_nSAF(self, plen='# AAs', **kwargs):
        '''
        '''

### check to see if user provided new data for cleaning
        if 'first' in kwargs.keys():
            if hasattr(self, 'first_data'):
                print('\nWarning: Copying over data provided in experiment() call.')
            self.first_data = kwargs['first']            
            
        if 'second' in kwargs.keys():
            if hasattr(self, 'second_data'):
                print('\nWarning: Copying over data provided in experiment() call.')
            self.second_data = kwargs['second']

### allow user to provide experimental data            
#         if 'experi' in kwargs.keys():
#             if hasattr(self, 'experi_data'):
#                 print('\nWarning: Copying over data provided in experiment() call.')
            
### check to make sure data provided matches `first` & `second` datasets
#             if kwargs['experi'].shape[0] != self.first_data.shape[0]:
#                 raise Exception('\nExperimental data provided does not match treatment data. Exiting...')

#             self.experi_data = kwargs['experi']

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
        
### do some quick plotting for diagnostic purposes
#         df_plot = self.experi_data[self.first_cols + self.second_cols].melt(var_name='Sample', value_name='Count')
#         g = sns.displot(x='Count', hue='Sample', data=df_plot, kind='kde', lw=2)
    
#         g.fig.set_size_inches(12,6)
#         sns.move_legend(obj=g, loc='upper right', bbox_to_anchor=(0.70, 1.0))
        
#         plt.show()
        
        return(self)
#----------------------------------------------------------------------------#
 

#----------------------------------------------------------------------------#
#     def plot_distr(self, xlabel=None, first=None, second=None):
    def plot_distr(self, log=True, testdata1=None, testdata2=None):

### define labels
        if log:
            xlabel = '$ln(NSAF)$'
            title = 'Log-Transformed Data Distributions'
        else:
            xlabel = 'NSAF'
            title = 'Normalized Data Distributions'

### define data for plotting        
        if not testdata1:
            if not hasattr(self, 'testdata1'):
                raise Exception('\nWarning: No statistical test performed. Cannot plot data. Exiting...')
            else:
                testdata1 = self.testdata1.copy()
                testdata1.columns = self.first_cols
        
        if not testdata2:
            if not hasattr(self, 'testdata2'):
                raise Exception('\nWarning: No statistical test performed. Cannot plot data. Exiting...')
            else:
                testdata2 = self.testdata2.copy()
                testdata2.columns = self.second_cols

        df_plot = testdata1.join(testdata2, how='outer')
        df_plot = df_plot.melt(var_name='Sample', value_name='Value')

### generate plot
        fig, ax = plt.subplots(1, 1, figsize=(12,8))

        sns.boxplot(x='Value', y='Sample', data=df_plot, whis=[0, 100], width=.6, palette='vlag', ax=ax)
        sns.stripplot(x='Value', y='Sample', data=df_plot, size=4, color='.3', linewidth=0, ax=ax)

        ax.set_xlabel(xlabel)
        ax.set_ylabel('')
        ax.set_title(title)
        plt.show()
        del(df_plot)
        
        return(fig)

#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### run statistical test for determining differential expression (t-test)
    def ttest(self, alpha=0.05, correction='BH', labels=None, log=True, **kwargs):
        '''

        Parameters
        ----------
        alpha : float, required
            Significance threshold for t-test, not adjusted for multiple comparisons. The default is 0.05.
        correction : string, optional... LIKELY A PROBLEM FOR THE VOLCANO PLOT
            Define correction for multiple comparison testing. The default is 'BH' (Benjamini-Hochberg, which right now is the only option).
        labels : string, optional
            Provides a list of labels to carry along for easier interpretation. The default is None.
        log : boolean, optional
            Choose whether to perform test on ln(NSAF) or NSAF. The default is True. 
        **kwargs : pd.DataFrame, optional
            Can provide new data for cleaning using 'first' & 'second' options, like in experiment() method. 

        Returns
        -------
        self.

        '''
        
### check to see if user provided new data for cleaning
        if 'first' in kwargs.keys():
            if hasattr(self, 'first_data'):
                print('\nWarning: Copying over data provided in experiment() call.')
            self.first_data = kwargs['first']            
            
        if 'second' in kwargs.keys():
            if hasattr(self, 'second_data'):
                print('\nWarning: Copying over data provided in experiment() call.')
            self.second_data = kwargs['second']
        
### assign attributes for future use
        self.alpha = alpha
        self.labels = labels
                
### get nSAF columns for statistical tests & plotting later
        self.l_nsafcols1 = [col for col in self.first_data.columns if 'nSAF' in col]
        self.l_nsafcols2 = [col for col in self.second_data.columns if 'nSAF' in col]
        
### define data for performing t-test... ln(NSAF) unless directed otherwise by user
        if log:
            self.testdata1 = np.log(self.first_data[self.l_nsafcols1])
            self.testdata2 = np.log(self.second_data[self.l_nsafcols2])
        else:
            self.testdata1 = self.first_data[self.l_nsafcols1]
            self.testdata2 = self.second_data[self.l_nsafcols2]

### run t-test
        ttest = stats.ttest_ind(a=self.testdata1, b=self.testdata2, alternative='two-sided', axis=1)
        l_pvals = ttest[1]

### plot test data distributions
#         self.plot_distr()

### join data for further analysis
        self.results = self.first_data.join(self.second_data)

### calculate fold change & ad to results frame
        self.results['fold-change'] = self.second_data[self.l_nsafcols2].mean(axis=1) / self.first_data[self.l_nsafcols1].mean(axis=1)

### add p-values
        self.results['p-value'] = l_pvals

### implement Benjamini-Hochberg procedure to correct for multiple comparisons... BREAK THIS OUT??
        if correction == 'BH':
            self.results['p-adj'] = smsm.multipletests(self.results['p-value'], alpha=alpha, method='fdr_bh')[1] # now this should be compared against alpha
            print(f'\nBonferroni-adjusted significance level:\t{smsm.multipletests(self.results["p-value"], alpha=alpha, method="fdr_bh")[3]:.5f}\n(just a ballpark check, Benjamini-Hochberg used here)')
        else:
            print('\nNo correction made for multiple comparisons...')

### add labels to frame if included            
        if labels:
            self.results = self.experi_data[[labels]].join(self.results, how='right')

        return(self)
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### function for calculating q-values
    def est_qval(self, results=None, pval='p-value'):
    
### raise issue if p-values haven't been generated yet
        if not results:
            if not hasattr(self, 'results'):
                raise Exception('\nData containing p-values required. Either provide using `results` parameter or run a statistical test first.\n\nExiting...')
            else:
                results = self.results

### new dataframe for converting p- to q-values
        df_q = results[[pval]].sort_values(pval).copy()
    
### define m as number of features under consideration
        m = df_q.shape[0]

### create a grid of hyperparameters, "lambda"
        h_grid = np.arange(0,0.96,0.01)

### estimate the expected proportion of null values, "pi0"
        l_pi = [(np.sum(df_q[pval] > h) / (m * (1 - h))) for h in h_grid]

### generate cubic spline of pi versus lambdas... NOT NATURAL, DESPITE OPTION; NOT SMOOTHED; POOR EXTRAPOLATION
#         spline_pi = scint.CubicSpline(x=h_grid, y=l_pi, bc_type='natural')


### train isotonic regression to estimate pi0... CURRENTLY UNTESTED WORKAROUND DUE TO LACK OF NORMAL CUBIC SPLINE IN PYTHON
        regress_pi = skliso.IsotonicRegression(increasing=False).fit(h_grid, l_pi) # WHAT ABOUT A NORMAL POLYNOMIAL REGRESSION VIA SKLEARN?
        pi0 = pd.Series(regress_pi.predict(h_grid)).value_counts().idxmax() # in essence, calulate the mode of the regression prediction
        print(f'\nExpected proportion of null features:\t\t{pi0:.2f}\n(used for calculating q-values)')


### quick diagnostic plot for estimation of pi0 required in calculating q-values... BREAK THIS OUT INTO SEPARATE METHOD?
        df_plot = pd.DataFrame.from_dict({'lambda':h_grid, 'pi0':l_pi, 'pi0-model':regress_pi.predict(h_grid)})
    
        fig, ax = plt.subplots(1, 1, figsize=(10,6))
        sns.scatterplot(x='lambda', y='pi0', data=df_plot, label='Data', ax=ax)
        sns.lineplot(x='lambda', y='pi0-model', data=df_plot, label='Model', ax=ax)
    
        ax.set_xlim([0,1])
        ax.set_ylim(top=1)
        ax.set_xlabel('Lambda')
        ax.set_ylabel('$\hat{\pi}_{0}$')
        ax.set_title('Estimation of the Proportion of Null Features ($\hat{\pi}_{0}$) using Isotonic Regression')
        plt.show()
        del(df_plot)


### calculate q-values
#         pm = pi0 * df_q.iloc[-1]['p-value']
        l_p = df_q['p-value'].tolist()[::-1]
        l_q = []
        
        for i,p in enumerate(l_p):
            if i == 0:
                l_q.append(pi0*m*p/(m-i)) # m-i instead of i b/c we reversed the order of p-values in l_p
            else:
                l_q.append(np.min([pi0*m*p/(m-i), l_q[i-1]])) # l_q[i-1] for comparison vice l_q[i+1] b/c we reversed the order of p-values in l_q
        
        df_q = df_q.sort_values('p-value', ascending=False)
        df_q['q-value'] = l_q
        df_q = df_q.sort_values('p-value')
        
        self.results['q-value'] = df_q['q-value']
        
        return(self)
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### q-value versus p-value
    def plot_pvq(self, results=None, pval='p-value', qval='q-value'):
 
### raise issue if q-values haven't been generated yet
        if not results:
            if not hasattr(self, 'results'):
                raise Exception('\nData containing p-values required. Either provide using `results` parameter or run a statistical test first.\n\nExiting...')
            else:
                df_plot = self.results
        else:
            df_plot = results

### generate plot
        df_diag = pd.DataFrame.from_dict({'x':np.arange(0,1.01,0.01), 'y':np.arange(0,1.01,0.01)})
    
        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        sns.lineplot(x=pval, y=qval, data=df_plot, label='q v. p', lw=4, ax=ax)
        sns.lineplot(x='x', y='y', data=df_diag, ls=':', label='y = x', ax=ax)
    
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        ax.legend(loc='lower right')
        ax.set_title('q-value v. p-value')
    
        return(fig)
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### # of significant genes versus q-value
    def plot_sigvq(self, results=None, pval='p-value', qval='q-value'):

### raise issue if q-values haven't been generated yet
        if not results:
            if not hasattr(self, 'results'):
                raise Exception('\nData containing p-values required. Either provide using `results` parameter or run a statistical test first.\n\nExiting...')
            else:
                df_q = self.results
        else:
            df_q = results

### generate plot
        l_qthresh = np.arange(0,0.11,0.0001)
        l_signif = [np.sum(df_q[qval] < qt) for qt in l_qthresh]
        df_plot = pd.DataFrame.from_dict({qval:l_qthresh, 'count':l_signif})
  
        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        sns.lineplot(x=qval, y='count', data=df_plot, lw=4, ax=ax)
    
        ax.set_xlim(np.min(l_qthresh), np.max(l_qthresh))
        ax.set_ylim(bottom=0)
        ax.set_ylabel('# of Significant Features')
        ax.set_title('Number of Signficant Features v. q-value')
        
        return(fig)
#----------------------------------------------------------------------------#


### run statistical test for determining differential expression (LIMMA)... GONNA BE TOUGH

#----------------------------------------------------------------------------#


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
            print('\nWorking on volcano plot labels. Can take a minute...')
            for idx in range(df_volc[df_volc['p-Sig'] == 1].shape[0])[:100]:
                x = 0.075+df_volc.loc[idx, 'log2(FC)']
                y = df_volc.loc[idx, '-log10(p)']
                s = df_volc.loc[idx, labels]
                plt.text(x=x, y=y, s=s, fontdict=dict(color='k', size=4))
        
        ax.set_xlim([xmin, xmax])
        ax.set_xlabel('$log_{2}$(Fold Enrichment)')
        ax.set_ylim([ymin, ymax])
        ax.set_ylabel('$-log_{10}$(p-value)')
#         plt.title(f'Volcano Plot ({name_outp})')
        plt.title(f'Volcano Plot')
        plt.show()
        
        return(fig) # return figure so user can save as s/he wants
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### define function for generating data formatted for over-representation analysis (ORA)
    def prep_ora(self, data=None, results=None, effect='fold-change', thresh=2.0):

### get reference data
        if data:
            df_ref = data.copy()
        else:
            df_ref = pd.DataFrame(self.experi_data[self.first_cols + self.second_cols].dropna(how='all').index.tolist(), columns=['Feature'])
    
### get results 
        if results:
            df_ora = results.sort_values(by=effect, ascending=False).copy()
        else:
            df_ora = self.results.sort_values(by=effect, ascending=False).copy()

### get over/under-represented features given threshold
        df_ora = pd.DataFrame(df_ora[(df_ora[effect] >= thresh) | (df_ora[effect] <= 1/thresh)].index.tolist(), columns=['Feature'])
        
        return(df_ora, df_ref)
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### define function for generating data formatted for set enrichment analysis (SEA)
    def prep_sea(self, effect='fold-change'):

### rank data by effect, lowest to highest
        df_sea = self.results.sort_values(by=effect)[[effect]].reset_index()
    
        return(df_sea)
#----------------------------------------------------------------------------#



##### CREATE CLASS FOR RUNNING STATISTICAL ENRICHMENT ANALYSIS
#----------------------------------------------------------------------------#
### create `enrichment` class for statistical enrichment analysis
class enrichment:
    '''
    Handles & analyzes proteomics data. 
    
    REQUIREMENTS:
        1. Mass spectrometry (MS) data w spectral counts. 
    '''
    
#----------------------------------------------------------------------------#
### initialization parameters... SHOULD I INCLUDE THESE PARAMETERS W THE EXPERIMENT FUNCTION INSTEAD?? THINK ABOUT HOW THIS WILL BE USED (NOT AS SINGLE SCRIPT, LIKE DONE BELOW)
    def __init__(self, path_data=None, accession=0):
        '''

        Parameters
        ----------
        path_geneExp : string, required
            Pathway to expression data, formatted w two columns representing protein accessions & fold-enrichment. The default is None.

        Returns
        -------
        None.

        '''
        self.geneExp = path_data
        try:
            self.seadata = pd.read_csv(self.geneExp, sep='\t', header=None) # GOING TO NEED TO ACCOUNT FOR OTHER DATA TYPES
        except:
            raise Exception('\nCannot open expression data. Exiting...')
            
        try:
            self.seagenes = self.seadata[accession].tolist()
        except:
            raise Exception('\nCannot find accession data. Exiting...')

### dictionary of pathways one can use... simplifies the process a bit for the user
        self.dict_annotset = {
            'GOMF':'GO:0003674', # GO molecular functions
            'GOBP':'GO:0008150', # GO biological process
            'Panther':'ANNOT_TYPE_ID_PANTHER_PATHWAY',
            'Reactome':'ANNOT_TYPE_ID_REACTOME_PATHWAY'
            }
        
### instantiate dictionary for storting protein names... meant to speed things up
        self.dict_protnames = {}
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### run set enrichment analysis via PANTHERDB
#     def run_sea(self, annotset='ANNOT_TYPE_ID_PANTHER_PATHWAY', organism=9606, correction='FDR', cutoff=None):
    def run_sea(self, annotset='Panther', organism=9606, correction='FDR', cutoff=None):
        '''
        
        Parameters
        ----------
        path_geneExp : string, required
            Path to .csv file containing enrichment data. As of now, must be a saved list of genes/proteins w two columns (accession number & fold-enrichment). The default is None.
        annotset : string, optional
            Defines which annotation set to run analysis against. The default is 'ANNOT_TYPE_ID_PANTHER_PATHWAY', representing the PANTHER pathway set.
        organism : integer, optional
            Identifies organism for analysis. The default is 9606, or human. 
        correction : string, optional
            Correction for multiple hypothesis testing. The default is 'FDR', which controls false discovery rate. 
        cutoff : float, optional
            User-defined cutoff for defining statistical enrichment. The default is None & will return all results.
    
        Returns
        -------
        None.
    
        '''
    
### the data you want to pull
        files = {'organism':organism, 'correction':correction, 'annotDataSet':self.dict_annotset[annotset], 'geneExp':open(self.geneExp, 'r')}
    
### url for PANTHER SEA API
        url = 'http://pantherdb.org/services/oai/pantherdb/enrich/statenrich'
    
### check connection to API... it can be finnicky sometimes
        try:
            requests.head(url, timeout=5) 
        except:
            raise Exception('Cannot reach API. Exiting...')
    
### access API... requires POST call since you have to run enrichment analysis, not just pull data
        with requests.post(url, files=files, timeout=5) as response: 
    
### extract enrichment results
            df_enri = pd.DataFrame.from_dict(response.json()['results']['result'])
            for idx, row in df_enri.iterrows():
                try:
                    df_enri.loc[idx, 'id'] = row['term']['id']
                    df_enri.loc[idx, 'label'] = row['term']['label']
                except:
                    print(f'\nNo pathway data for index {idx}. Moving on...')
            
### format the data
            l_seacols = ['id', 'label', 'number_in_list', 'fdr', 'pValue', 'plus_minus']
            dict_seacols = {'id':'Pathway', 'label':'Description', 'number_in_list':'Count', 'fdr':'FDR', 'pValue':'p-value', 'plus_minus':'+/-'}
            df_enri = df_enri[l_seacols].rename(columns=dict_seacols).sort_values(by='FDR').copy()
                        
### get list of pathway IDs for mapping to proteins
            if cutoff:
                df_enri = df_enri[df_enri['FDR'] < cutoff].copy()

            print(f'\n{df_enri}\n')
            
### set new attributes for use later
            self.enrichedpaths = df_enri['Pathway'].tolist()                
            
        return(df_enri)
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
### match proteins in sea data to enriched pathways & pull down related data from PANTHERDB
#     def match_sea(self, annotset='ANNOT_TYPE_ID_PANTHER_PATHWAY', organism=9606, df_sea=None, l_seagenes=None, l_enripaths=None):
    def match_sea(self, annotset='Panther', organism=9606, df_sea=None, l_seagenes=None, l_enripaths=None):

### use enrichment data resulting from run_sea() method, if not provided
        if not df_sea:
            df_sea = self.seadata

### use enrichment data resulting from run_sea() method, if not provided
        if not l_seagenes:
            l_seagenes = self.seagenes
            
### use enrichment data resulting from run_sea() method, if not provided
        if not l_enripaths:
            l_enripaths = self.enrichedpaths

### pulls gene data for those in enriched set... MAX 1000 @ A TIME, WILL NEED TO ACCOUNT!!!
        genes = '%2C'.join(l_seagenes) # %2C serves as the delimiter for url
        url = 'http://pantherdb.org/services/oai/pantherdb/geneinfo?' + f'geneInputList={genes}&organism={organism}'
    
### get list of enriched pathways
#         l_enripaths = df_enri['Pathway'].tolist()
    
### access API... only pulling data for proteins in enrichment set
        with requests.get(url) as response:
    
### count for determining fraction of instances mapped to particular annotation set
            i_set = 0
            i_tot = 0
        
### start dictionary for mapping proteins to pathways
            dict_pp = {}
        
### filter data to proteins in enrichment list
            for gene in response.json()['search']['mapped_genes']['gene']:
        
### check to ensure data completeness...
                if ('annotation_type_list') not in gene.keys():
                    print(f'No annotations on {gene["accession"]}. Moving on...')
                    pass
    
### if available, check annotations    
                else:    
                    if type(gene['annotation_type_list']['annotation_data_type']) is not list: # in cases where only one entry, not a list
                        l_gene = [gene['annotation_type_list']['annotation_data_type']]
                    else:
                        l_gene = gene['annotation_type_list']['annotation_data_type']
        
### check for match w annotation data set (e.g., PANTHER Pathways or GO Biological Process) 
                    for annot in l_gene:
                        if annot['content'] == self.dict_annotset[annotset]:
                            
                            if type(annot['annotation_list']['annotation']) is not list:
                                l_annot = [annot['annotation_list']['annotation']]
                            else:
                                l_annot = annot['annotation_list']['annotation']
    
### add data (right now just UniProtID) to pathway dictionary if so                    
                            for a in l_annot:                        
                                if a['id'] in l_enripaths:
                                    uniprotid = gene['accession'].split('|')[-1].split('=')[-1]
                            
                                    if a['id'] not in dict_pp.keys():
                                        dict_pp[a['id']] = [uniprotid]
                                    else:
                                        dict_pp[a['id']].append(uniprotid)
    
### update count for number of proteins in set w matched pathways                            
                            i_set += 1
    
### update count for total number of proteins queried            
                i_tot += 1
    
### print statistic for user... tells you fraction of proteins matched to a pathway
            print(f'\n{100*(1-(i_set/i_tot)):.1f}% protein IDs not found...')
            
### convert dictionary to long form dataframe    
        df_protid = pd.DataFrame.from_dict(dict_pp, orient='index').melt(ignore_index=False).reset_index().dropna()[['index', 'value']].rename(columns={'index':'Pathway', 'value':'Accession'}).sort_values(by=['Pathway', 'Accession'])
    
        return(df_protid)



### define function for pulling protein name & information from UniProt API
    def get_protname(self, accession=None, organism_id=9606):
        '''
    
        Parameters
        ----------
        accession : string, required
            The UniProt accession number for the protein of interest. The default is 'P60709'.
        organism_id : numeric, optional
            The organism ID for searching proteins. The default is 9606, the ID for humans.
    
        Returns
        -------
        Right now, only the protein name. 
    
        '''
    
### ensure accession provided before anything else
        if not accession:
            raise Exception('\nWarning: No protein accession provided. Exiting...')
            
### pull protein name from local dictionary if already searched for...
        if accession in self.dict_protnames.keys():
            message = f'{accession} already identified'
            protname = self.dict_protnames[accession]
            
### ...otherwise, pull from PANTHER
        else:
            message = f'Working on {accession}'
        
### define url for pulling protein name from API... IS THERE AN EASIER WAY? 
            url = f'https://rest.uniprot.org/uniprotkb/search?query=accession:{accession}&organism_id:{organism_id}&columns=protein_name'    
    
### request uniprot API & print status
            with requests.get(url) as response:
                if response.status_code != 200:
                    raise Exception(f'\nCaution: Link to UniProt API failed for {accession}. Exiting...')
            
### pull data of interest... right now, just uniprotid (for confirmation) & protein name
                protname = response.json()['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
                self.dict_protnames[accession] = protname
        
### message for user
        if len(protname) > 50:
            protname_print = protname.split(',')[0]
        else:
            protname_print = protname
            
        print(f'{message}...{protname_print}')
            
### close request... not sure this is necessary but can't hurt
#             response.close()
        
        return(protname)
















