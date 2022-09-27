# Omics
Module for analyzing omics data written in Python.
Will provide classes & functions for proteomics analysis, to start. 

## `omics.py`
The beginnings of a module designed to provide functions for interpreting & analyzing omics data. 
**Currently houses the `dpea` class which enables differential protein expression analyses using mass spectrometry (MS) data.** The `dpea` class requires your MS data as initial input. 
The following methods are included:
- `.experiment()` defines experimental parameters & formats data appropriately. This is the preferred starting point for any `dpea` workflow, as it makes running the rest of the methods simpler. Methods are being developed to enable independent & modular usage, however. 
- `.clean()` cleans up your MS data. Deals with missing data according to your preference & adds pseudocounts to ensure no calculations wind up as infinite. 
- `.norm_nSAF()` calculates normalized spectral abundance factors from spectral counts, thus normalizing your data. The calculation is as defined by [this 2012 paper by McIlwain et al.](https://doi.org/10.1186/1471-2105-13-308)
- `.ttest()` performs a two-sided t-test for each feature provided across samples. Adjustments are made for multiple comparisons using the `correction` parameter, which as of now only includes the Bengamini-Hochberg method (but will be expanded!). **p-values & adjusted p-values are provided in results.**
- `.est_qval()` calculates the q-values, or the false discovery rates, according to [this 2003 paper by Storey & Tibshirani](https://doi.org/10.1073/pnas.1530509100). 
- `.plot_distr()` ...
- `.plot_pvq()` ...
- `.plot_sigvq()` ...
- `.volcano()` ...
- `.prep_ora()` prepares your data for representation analyis using [PANTHER](http://pantherdb.org/). Yields two data frames: one containing the proteins which meet a user-defined threshold for differential expression & another which contains all the proteins detected in your samples. 
- `.prep_sea()` prepares your data for enrichment analyis using [PANTHER](http://pantherdb.org/). Yields a single data frame ranked by differential expression. 

The `dpea` class requires your proteomics dataset along w some parameters in calling. 
The primary method is called `experiment`, which will transform your data & perform statistical tests according to defaults. 
Methods are currently being reworked to enable more independence. 
Highly suggest using `params.json` file for input parameters (see `example.json` provided); doing so provides reproducibility & is generally simpler. 

Working on a new `enrichment` class which builds on `dpea`. Goal is for it to implement via [PANTHER API](http://pantherdb.org/services/openAPISpec.jsp).
The following methods are in progress:
- `run_sea()` ...
- `analyze()` ... should probably rename this to something more specific.

## Exp
Contains json files w script parameters for each experiment. 

## Env
Contains conda environment yml file.  

## Params
Contains two sample parameter files. Not required, but highly recommended! Simplifies the analysis & improves reproducibility. 