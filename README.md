# Omics
Module for analyzing omics data written in Python.
Will provide classes & functions for proteomics analysis, to start. 

## `omics.py`
The beginning of a module designed to provide functions for interpreting & analyzing omics data. 
Currently holds only the `dpea` class, containing methods for running differential protein expression analyses. 
The `dpea` class requires your proteomics dataset along w some parameters in calling. 
The primary method is called `experiment`, which will transform your data & perform statistical tests according to defaults. 
Methods are currently being reworked to enable more independence. 
Highly suggest using `params.json` file for input parameters (see `example.json` provided); doing so provides reproducibility & is generally simpler. 


## Exp
Contains json files w script parameters for each experiment. 

## Env
Contains conda environment yml file.  
