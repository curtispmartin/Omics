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

### set analysis parameters
alpha = float(dict_params['Parameters']['alpha']) # confidence level
fcthresh = float(dict_params['Parameters']['fcthresh']) # threshold for relevant fold-change
plen = dict_params['Metadata']['Lengths']


### instantiate class object using data
first = dict_params['Treatments']['1']
second = dict_params['Treatments']['2']
names = 'Accession'
test = omics.dpea(df=df)
test = test.experiment(first=first, second=second, names=names)

df_first = test.first_data
df_second = test.second_data
df_results = test.results

fig = test.volcano(labels='GenSymbol')
# fig.savefig(f'{name_outp}.png', bbox_inches='tight', dpi=300)

    

