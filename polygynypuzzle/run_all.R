# This code will replicate all tables and figures from the paper
# Be sure to run each line in order
 rm(list=ls())

# Set your own path, everything else is relative
 path1 <- 'C:\\Users\\Mind Is Moving\\Dropbox\\Completed and Published Projects\\1-papers\\31-Polygyny Puzzle\\Interface\\Revisions\\SupplementaryMaterials'
 setwd(path1)

source('./Code/describe_data.R')
source('./Code/fit_models.R')
source('./Code/explore_model.R')
source('./Code/extrapolate_fits.R')
source('./Code/explore_extrapolations.R')
source('./Code/explore_extrapolations_with_extra_populations.R')
source('./Code/check_model.R')




