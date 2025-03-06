library(tidyverse)
library(readxl)

omnipath_network <- read_csv('../data/OmniPath/2023_10_11_KinaseDataOmniPath.csv')
omnipath_interactions <- unique(paste0(omnipath_network$enzyme, '_', omnipath_network$TARGET_UP_ID))
# Load the data
baseline_KIN <- read_tsv('../results/baselineKIN_comparisons/baseline_KIN/buljan2020_comparison_KIN.tsv')
baseline_KIN <- baseline_KIN %>% filter(Edge_type == 'KS')

data <- read_excel('../data/competitors/buljan2020_interactions.xlsx', na = c('NA', ''))
colnames(data) <- data[2, ]
data <- data %>%
  filter(!row_number() %in% c(1, 2)) %>%
  filter(!is.na(wdscore), !is.na(gfpratio)) %>%
  mutate(wdscore = as.numeric(wdscore), gfpratio = as.numeric(gfpratio)) %>%
  filter(wdscore >= 73.6, gfpratio >= 18.4)


# filtering for theoretically recoverable interactions
baseline_KIN <- baseline_KIN %>% filter(Source %in% data$Bait_id)
data <- data %>% filter(Bait_id %in% baseline_KIN$Source)

# filtering for top 15 interactions
baseline_KIN <- baseline_KIN %>% group_by(Target) %>% slice_max(n = 15, order_by=-percentileRank) %>% ungroup()

buljan_interactions <- c(unique(paste0(data$Bait_id, '_', data$Protein_id)))
baseline_KIN_interactions <- unique(paste0(baseline_KIN$Source, '_', baseline_KIN$Target_Uniprot))

print(length(buljan_interactions))
print(length(baseline_KIN_interactions))
method_intersection <- length(intersect(buljan_interactions, baseline_KIN_interactions))
baseline_op_intersection <- length(intersect(baseline_KIN_interactions, omnipath_interactions))
buljan_op_intersection <- length(intersect(buljan_interactions, omnipath_interactions))
buljan_baseline_op_intersection <- length(intersect(intersect(buljan_interactions, omnipath_interactions), baseline_KIN_interactions))

print(method_intersection)
print(baseline_op_intersection)
print(buljan_op_intersection)
print(buljan_baseline_op_intersection)

print(method_intersection / length(buljan_interactions))
print(baseline_op_intersection / length(baseline_KIN_interactions))
print(buljan_op_intersection / length(buljan_interactions))