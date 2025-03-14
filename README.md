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

## R
This folder contains R code to tidy the environmental and eDNA data, fit models, and plot preliminary results. Files within a folder are prefixed with a number to indicate which order they should be run in, if there is no prefix, files can be run in any order. 

#### get_envfeatures
Contains files to retrieve predictive features from raw environmental data.

#### ebv_fit
Contains files to derive EBVs, generate model fitting scripts, and plot results of model fitting

#### Dardel 
Contains files to fit models on the dardel cluster

#### Misc
Miscellaneous files useful for the FinBio project

## Results
Results files from model fits and predictions

## Plots
Plotted figures.
