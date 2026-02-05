# finbio-ebv
## Introduction
This repository contains code to run analyses aimed at evaluating the use of large scale eDNA inventories in conjunction with earth observation (EO), to predict essential biodiversity variables (EBV's) for use in biodiversity impact reporting. The project aims to asses whether summary measures of biodiversity can be accurately predicted using common EO products and machine learning algorithms. The code derives several EBV's that summarise key components of invertebrate diversity using data from the Insect Biome Atlas project, collected across Sweden and Madagascar. Gradient boosted regression tree's (XGboost), are used to learn associations between these summary measures and environmental features derived from the ESA's copernicus landcover data, and ERA5 climate data. Cross-validation is performed to assess predictive performance within sites and at new sites. 

## Data
This repository contains only the tidied data files, raw eDNA and environmental data should be downloaded from [10.6084/m9.figshare.28590113](https://figshare.com/articles/dataset/Raw_data_for_FinBio_EBV_project/28590113), and copied into the 'data/' directory.
The folder contains three subfolders:
#### raw_envdata
Containing raw land cover (copernicus) and era5 measurements. 
#### IBA_data
Containing the files required to derive EBV's from the IBA data. 
#### tidy_data 
Containing the tidied environmental data and derived EBV's

Also included is a phylogenetic tree (species_level_tree.nwk), which is used for the derivation of a measure of phylogenetic diversity. 

'full_M.rds' is a list containing each EBV alongside the environmental features required for model fitting - this is the only data file required to fit the models, the other components are provided for clarity. 

## R
This folder contains R code to tidy the environmental and eDNA data, fit models, and plot preliminary results. Files within a folder are prefixed with a number to indicate which order they should be run in, if there is no prefix, files can be run in any order. 

#### get_envfeatures
Contains files to retrieve predictive features from raw environmental data.

#### ebv_fit
Contains files to derive EBVs from eDNA and trait data, generate model fitting scripts, and plot results of model fitting:

1_tidy_IBA_data.R - A script to do some initial tidying to the IBA community and meta data.  

2_derive_community_EBVs.R - A script to derive community level EBVs (Species richness and LCBD).  

3_derive_functional_EBVs.R - A script to derive functional eveness and dispersion.  

4_derive_genetic_EBVs.R - A script to derive genetic and phylogenetic diversity.  

5_build_ft_matrix.R - A script to build matrices of features and prediction targets for each EBV.  

6_generate_dardel_scripts.R - A script to programatically build model fitting scripts and shell scripts for jobs on the HPC cluster.  

7_dardel - A subfolder that contains data, .sh files and .R files to fit models.  

8_plot_loss.R - A script to plot results of cross-validation exercise.  

9_plot_ebv_trends.R - A script to visualise estimates of broad spatiotemporal trends in EBV's.

10_plot_rank_changes.R - A script to plot changes in assigned priority tiers to simulate a simple process for site prioritisation.


#### Misc
Miscellaneous files useful for the FinBio project

## Results
Results files from model fits and predictions

## Corporate report analysis
Excel document with the full list of reports reviewed for the corporate report analysis in the study.

## Plots
Plotted figures

# ------------------------------------ #
## Dardel setup
These are some basic instructions to set up this code to run on the dardel HPC cluster: https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system/dardel-1.1043529 

1) Submit a job to install R packages - these would not install in an interactive session so I installed via a batch job. An R script with the following should install everything required:


install.packages("tidyverse",repos='http://cran.us.r-project.org')
install.packages("xgboost",repos='http://cran.us.r-project.org')


2) Create an directory with two subfolders
   -'results'
   -'logfiles'
The model fitting code will store all results in the 'results' subfolder, and all output from dardel in the logfiles subfolder.

3) Three files are required to run a model -
   a) full_M.rds - the data file containing all EBV and environmental data. Each list component is an EBV (in column 2), and a series of environmental features.  
   b) The model file for a specific EBV * country combination, E.g. fit_model_madagascar_FD.R , is a model to run out CV exercise for functional dispersion (FD) metrics in madagascar.
   c) The batch file for a specific EBV * country combination - submit_model_madagascar_FD.sh. 






