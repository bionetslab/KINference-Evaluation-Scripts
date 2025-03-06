library(tidyverse)
library(readxl)
library(data.table)

omnipath_network <- read_csv('../data/OmniPath/2023_10_11_KinaseDataOmniPath.csv')
omnipath_interactions <- unique(paste0(omnipath_network$enzyme, '_', omnipath_network$TARGET_UP_ID))
# Load the data
baseline_KIN <- read_tsv('../results/Hyperparameters_old/Bouhaddou2023/baseline_KIN/Bouhaddou2023_n15_alpha0.9_KIN.tsv')
baseline_KIN <- baseline_KIN %>% filter(Edge_type == 'KS')

grnboost2_x0 <- read_tsv('../results/GRNBoost2/x0_network.tsv', col_names = F)
grnboost2_x1 <- read_tsv('../results/GRNBoost2/x1_network.tsv', col_names = F)
# Filtering for top 100,000 edges
# grnboost2_x0 <- grnboost2_x0[1:100000,] 
# grnboost2_x1 <- grnboost2_x1[1:100000,]

grnboost2_x0$Source <- sapply(strsplit(grnboost2_x0$X1, '_'), function(x){ x[1] })
grnboost2_x1$Source <- sapply(strsplit(grnboost2_x1$X1, '_'), function(x){ x[1] })
grnboost2_x0$Target <- sapply(strsplit(grnboost2_x0$X2, '_'), function(x){ x[1] })
grnboost2_x1$Target <- sapply(strsplit(grnboost2_x1$X2, '_'), function(x){ x[1] })
grnboost2_interactions <- unique(c(paste0(grnboost2_x0$Source, '_', grnboost2_x0$Target), paste0(grnboost2_x1$Source, '_', grnboost2_x1$Target)))

baseline_KIN_interactions <- unique(paste0(baseline_KIN$Source, '_', baseline_KIN$Target_Uniprot))

print(length(grnboost2_interactions))
print(length(baseline_KIN_interactions))
method_intersection <- length(intersect(grnboost2_interactions, baseline_KIN_interactions))
baseline_op_intersection <- length(intersect(baseline_KIN_interactions, omnipath_interactions))
grnboost2_op_intersection <- length(intersect(grnboost2_interactions, omnipath_interactions))

print(method_intersection)
print(baseline_op_intersection)
print(grnboost2_op_intersection)
print(length(intersect(intersect(grnboost2_interactions, omnipath_interactions), baseline_KIN_interactions)))

print(baseline_op_intersection / length(baseline_KIN_interactions))
print(grnboost2_op_intersection / length(grnboost2_interactions))