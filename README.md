# Omics
Module for analyzing omics data written in Python.
Will provide classes & functions for proteomics analysis, to start. 

## `omics.py`
Currently performs differential protein expression analysis using mass spectrometry data. 
Highly suggest using `params.json` file as input; doing so provides reproducibility & is generally simpler. 
However code can be used in command line w the following options:
- `-i` indicates a json input data file.
- `-f` points to the mass spec data file & indicates further user options required. 
Working to generalize &, eventually, convert to python module. 

## Exp
Contains json files w script parameters for each experiment. 

## Env
Contains conda environment yml file.  
