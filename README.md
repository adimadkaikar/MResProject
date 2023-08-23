# MResProject

### This is a repository containing all the code and data used for the project 'Predicting the thermal niche of a ubiquitous bacterium using whole genome sequence'

This repository contains the code, and data required to recreate the simulations perfomed to quantify the growth rate of E. coli across a range of temperatures. It contains three directories: 'code', 'data', 'results', 'models'
## Languages 
* R(v4.1.2)
* Python(v3.10.0)

## Dependencies 
* COBRApy(v0.26.3)
* Gurobi LP solver(v9.5.1rc2)
* Scikit-learn
* tidyverse
* cowplot

## Code

This directory contains all the code files used:
 
* GEMS.py -> For laoding the infput files in correct data strucutres
* LoopySimulation.ipynb -> Prediction of growth rates over all empirical using final parameter values
* onlyKcat.ipynb -> Addition of only enzyme activity information to ecGEM
* onlyNGAM -> Addition of only cell maintenence information to ecGEM
* onlyTm -> Additon of intracellular enzyme concentration information to ecGEM
* plotting.R -> Creating plots required with simuaiton results
* plotting2.R -> Code for plotting the effet of each temperature variant factor
* plottingMedia.R -> Creating plots for changing media simulations
* plottingInset.R -> Creating inset plots
* relaxingBounds.ipynb -> Temperature dependent simulation by relaxing bounds on reactions denoting enzyme abundance
* simulateGEM.ipynb -> Temperature dependent simulation with initial paramters
* SimulationLB.ipynb -> Temperature dependent simulation with final parameters and LB medium
* SimulationMM.ipynb -> Temperature dependent simulation with final parameters and Minimal medium
* TmToptCor.R -> Code for checking the correlation in Tm and Topt in initial parameters
* ToptFixed.ipynb -> Temperature dependent simulation by increasing the values for Topt
* ToptTmTopt.ipynb -> Multiple temperature dependent simulations following the algorithm presented
* TPCPlot.R -> Code to plot TPC in Fig. 1b
* wigCp.ipynb -> Code to check the effect of heat capacity change on model predictions
* wigTm.ipynb -> Code to check the effect of enzyme melting temperaure on model predictions
* wigCp.ipynb -> Code to check the effect of temperautre of maximum enzyme activity on model predictions
* etcpy/etc.py -> This file contains the temperature dependence modification code
* etcpy/tempDepCode.py -> This file contains the code for parameter sampling

## Data

This directory contains the data required for simulations:

* BestParamsTopt.csv -> Best parameter estimates achieved by increasing Topt
* BestParamsTopt129.csv -> Best parameter estimate
* ExpGrwoth.csv -> Empirical growth data used for model prediction and validation
* model_enzyme_params_new_tageed.csv -> Initial parameter estimates
* model_enzyme_params_new.csv -> Initial parameter estimates
* model_enzyme_params.csv -> Initial parameter estimates
* media.csv -> File containing media information
## Models

This directory contains the ecGEM used for incorporating temperature dependence of enzyme function

## Results

This directory contains results of the simulations