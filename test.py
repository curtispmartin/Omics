#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 22:43:22 2022

@author: martincup

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
# import scipy.interpolate as scint # for fitting natural cubic spline in calculating q-values
import sklearn.isotonic as skliso # utilize the isotonic regression function in place of natural cubic spline which isn't available
import statsmodels.stats.multitest as smsm # this could be a challenge... worth it? 

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/Users/martincup/Research/Omics') # add Omics directory to python path
import omics # work in progress!


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
path_outp = os.path.join(os.getcwd(), 'Output', name_inpu)

### create output directory if none exists
if not os.path.exists(path_outp):
    os.makedirs(path_outp)


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
alpha = float(dict_params['Parameters']['alpha']) # confidence level
exp = exp.ttest(alpha=alpha, correction='BH', labels='GenSymbol')

### try running all steps in compound method
# exp = omics.dpea(df=df).experiment(first=first, second=second, names=names).clean(nfloor=10, pseudocount=1).norm_nSAF(plen='# AAs').ttest(alpha=0.05, correction='BH', labels='GenSymbol')

### generate volcano plot
fcthresh = float(dict_params['Parameters']['fcthresh']) # threshold for relevant fold-change
fig = exp.volcano(fcthresh=fcthresh, labels='GenSymbol')
fig.savefig(os.path.join(path_outp, f'{name_outp}-volcano.png'), bbox_inches='tight', dpi=300)

### generate q-values & some diagnostic plots
exp = exp.gen_qval(pval='p-value')

fig = exp.plot_pvq()
fig.savefig(os.path.join(path_outp, f'{name_outp}-pvq.png'), bbox_inches='tight', dpi=300)

fig = exp.plot_sigvq()
fig.savefig(os.path.join(path_outp, f'{name_outp}-sigvq.png'), bbox_inches='tight', dpi=300)

### prepare data for ORA
df_ora, df_ref = exp.prep_ora()
df_ora.to_csv(os.path.join(path_outp, f'{name_outp}-ora.txt'), index=False, header=False)
df_ref.to_csv(os.path.join(path_outp, f'{name_outp}-ref.txt'), index=False, header=False)

### prepare data for SEA
df_sea = exp.prep_sea()
df_sea.to_csv(os.path.join(path_outp, f'{name_outp}-sea.txt'), sep='\t', index=False, header=False)

### save data for fidelity
df_outp = exp.results.copy()
df_outp.loc[df_outp[df_outp['q-value'] < alpha].index, 'significance'] = 1
df_outp['significance'] = df_outp['significance'].fillna(0)
df_outp = df_outp.sort_values(by=['significance', 'fold-change'], ascending=False)
df_outp.to_csv(os.path.join(path_outp, f'{name_outp}-processed.csv'))







    

