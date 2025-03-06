library(tidyverse)
library(readxl)
library(data.table)

omnipath_network <- read_csv('../data/OmniPath/2023_10_11_KinaseDataOmniPath.csv')
omnipath_interactions <- unique(paste0(omnipath_network$enzyme, '_', omnipath_network$TARGET_UP_ID))
# Load the data
baseline_KIN <- read_tsv('../results/baselineKIN_comparisons/baseline_KIN/invergo2020_comparison_KIN.tsv')
baseline_KIN <- baseline_KIN %>% filter(Edge_type == 'KS')

data <- read_excel('../data/competitors/invergo2020_interactions.xlsx', na = c('NA', ''))
colnames(data) <- data[1, ]
data <- data %>% filter(!row_number() %in% c(1)) %>%
    filter(!is.na(Mean_posterior_probability_of_regulation)) %>%
    mutate(Mean_posterior_probability_of_regulation = as.numeric(Mean_posterior_probability_of_regulation)) %>%
    filter(Mean_posterior_probability_of_regulation >= 0.5)

# Load Biomart data
biomart <- read_table('../data/biomart/mart_export_human_uniprot.txt') %>%
  select(ID_1, version) %>%
  distinct() %>%
  rename(Gene_ID = ID_1, Uniprot_ID = version) %>%
  filter(!is.na(Uniprot_ID))

data$Transducing_kinase <- sapply(data$Transducing_kinase, function(x){ (biomart %>% filter(Gene_ID == x) %>% select(Uniprot_ID) %>% pull())[1] })
data$Substrate_kinase <- sapply(data$Substrate_kinase, function(x){ (biomart %>% filter(Gene_ID == x) %>% select(Uniprot_ID) %>% pull())[1] })

# filtering for theoretically recoverable interactions
baseline_KIN <- baseline_KIN %>% filter(Source %in% data$Transducing_kinase)
data <- data %>% filter(Transducing_kinase %in% baseline_KIN$Source)

# filtering for top 15 interactions
baseline_KIN <- baseline_KIN %>% group_by(Target) %>% slice_max(n = 15, order_by=-percentileRank) %>% ungroup()

invergo_interactions <- c(unique(paste0(data$Transducing_kinase, '_', data$Substrate_kinase)))
baseline_KIN_interactions <- unique(paste0(baseline_KIN$Source, '_', baseline_KIN$Target_Uniprot))

print(length(invergo_interactions))
print(length(baseline_KIN_interactions))
method_intersection <- length(intersect(invergo_interactions, baseline_KIN_interactions))
baseline_op_intersection <- length(intersect(baseline_KIN_interactions, omnipath_interactions))
invergo_op_intersection <- length(intersect(invergo_interactions, omnipath_interactions))
invergo_baseline_op_intersection <- length(intersect(intersect(invergo_interactions, omnipath_interactions), baseline_KIN_interactions))

print(method_intersection)
print(baseline_op_intersection)
print(invergo_op_intersection)
print(invergo_baseline_op_intersection)

print(method_intersection / length(invergo_interactions))
print(baseline_op_intersection / length(baseline_KIN_interactions))
print(invergo_op_intersection / length(invergo_interactions))