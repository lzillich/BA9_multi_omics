# Author: Eric Zillich, last change: 2024-04-04
# MOFA with adapted code from https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html

library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(psych)
library(ggpubr)

DF <- data.frame

setwd(dir0)

#Specify the python version that will be used for training the model
reticulate::use_python("/anaconda3/bin/python", required=TRUE)

# import input data
BA9 <- get(load("/MOFA_dir/input/BA9_meth_expr_Rdata"))

# Add names and bring into the same order
names(BA9) <- c("Meth","Expr")
BA9$Meth <- BA9$Meth[,order(as.numeric(colnames(BA9$Meth)),decreasing = F)]
BA9$Expr <- BA9$Expr[,order(as.numeric(colnames(BA9$Expr)),decreasing = F)]

  # Create MOFA object
  MOFAobject <- create_mofa(BA9)
  
  # Specify MOFA options
  data_opts <- get_default_data_options(MOFAobject)

  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- 10
  
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- 42
  train_opts$maxiter <- 10000

 
  # Train the model
  MOFAobject <- prepare_mofa(MOFAobject,
                             data_options = data_opts,
                             model_options = model_opts,
                             training_options = train_opts
  )
  
  MOFAobject <- run_mofa(MOFAobject, outfile="/MOFA_dir/output/MOFA_BA9_trained_model_new.hdf5",use_basilisk = F)
  saveRDS(MOFAobject,"/MOFA_dir/output/MOFA_BA9_trained_model_20240102.rds")



