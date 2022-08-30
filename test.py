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
import statsmodels.stats.multitest as smsm # this could be a challenge... worth it? 

import matplotlib.pyplot as plt
import seaborn as sns

import omics # work in progress!


### pull data for testing
path_data = 'Data/KB-Appendix17-P798-16-IPs-Summary-Comparison.xlsx'
df = pd.read_excel(io=path_data, sheet_name='Summary_Comparison')

path_params = 'Data/Params/20220815_22Rv1_UN_IACS.json'
dict_params = json.load(open(path_params))

### get IO variables
name_inpu = os.path.split(dict_params['Metadata']['Data'])[-1].split('.')[0]
name_sheet = dict_params['Metadata']['Sheet']
path_data = os.path.split(dict_params['Metadata']['Data'])[0]
path_inpu = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0], dict_params['Metadata']['Data'])
name_outp = dict_params['Metadata']['Name']


### create dpea object
test = omics.dpea(df=df)

### define experiment
first = dict_params['Treatments']['1']
second = dict_params['Treatments']['2']
names = 'Accession'
test = test.experiment(first=first, second=second, names=names)

### clean data
test = test.clean(nfloor=10, pseudocount=1)

### normalize data
plen = dict_params['Metadata']['Lengths']
test = test.norm_nSAF(plen=plen)

### run statistical test
alpha = float(dict_params['Parameters']['alpha']) # confidence level
test = test.ttest(alpha=alpha, correction='BH', labels='GenSymbol')

### try running all steps in compound method
# test = omics.dpea(df=df).experiment(first=first, second=second, names=names).clean(nfloor=10, pseudocount=1).norm_nSAF(plen='# AAs').ttest(alpha=0.05, correction='BH', labels='GenSymbol')

### generate volcano plot
fcthresh = float(dict_params['Parameters']['fcthresh']) # threshold for relevant fold-change
fig = test.volcano(fcthresh=fcthresh, labels='GenSymbol')
# fig.savefig(f'{name_outp}.png', bbox_inches='tight', dpi=300)





    

