

setwd("~/projects/lipid_converting_enzymes/")

library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
library(readr)

data_exp=read.csv('signatures_lm_xpr_pert_pc1removed.csv',sep = ',',header = T,row.names = 1, check.names = FALSE)
data_exp_t <- t(data_exp)


data(dorothea_hs, package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

tf_activities <- run_viper(data_exp, regulons, 
                           options =  list(method = 'none', minsize = 4, 
                                           eset.filter = FALSE, cores = 4, 
                                           verbose = FALSE))
tf_activities <- t(tf_activities)

write.csv(tf_activities,'dorothea_ABC_signatures_lm_xpr_pert_pc1removed.csv')