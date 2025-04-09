#' Simulation dependencies
#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 24JAN2025

# Libraries
library(tidyverse) # dataframe manipulation
library(nlme)      # linear mixed models
library(EnvStats)  # truncated log normal distribution
library(psych)     # tetrachoric correlation 
library(MASS)      # multivariate normal

# Scripts
source('src/help.R')            # helper functions
source('src/genResponse.R')     # fxn to generate responses given scenario
source('src/simScenarios.R')    # fxn to simulate given scenarios

