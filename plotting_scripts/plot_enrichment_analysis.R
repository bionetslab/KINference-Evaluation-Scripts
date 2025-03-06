library(tidyverse)
library(gprofiler2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

colors <- display.brewer.all(n = 3, colorblindFriendly = TRUE)

colors <- brewer.pal(n = 3, "Dark2")
# Bouhaddou 2023
bouhaddou2023_baseline_KIN <- read_tsv('../results/Hyperparameters/Bouhaddou2023/baseline_KIN/Bouhaddou2023_n15_alpha0.9_KIN.tsv')
bouhaddou2023_filtered_KIN <- read_tsv('../results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma1_beta0.4_PCST_filtered_KIN.tsv')

data <- read_tsv('../data/Bouhaddou2023/intensities_Mock10h.tsv')
bg <- unique(sapply(strsplit(data$Protein, '_'), function(x) { x[1] }))

kinases <- c(
  read_tsv('../data/kinase_data/serine_threonine_kinase_name_mappings.tsv')[['ACC#']],
  read_tsv('../data/kinase_data/tyrosine_kinase_name_mappings.tsv')[['ACC#']]
)

query_baseline_KIN <- unique(bouhaddou2023_baseline_KIN[which((bouhaddou2023_baseline_KIN$Edge_type == 'KS') & !(bouhaddou2023_baseline_KIN$Target_Uniprot %in% kinases)),]$Target_Uniprot)
bouhaddou2023_gostres_baseline_KIN <- gost(
  query = query_baseline_KIN,
  organism = 'hsapiens',
  sources = c('KEGG')
)

query_filtered_KIN <- unique(bouhaddou2023_filtered_KIN[which((bouhaddou2023_filtered_KIN$Edge_type == 'KS') & !(bouhaddou2023_filtered_KIN$Target_Uniprot %in% kinases)),]$Target_Uniprot)
bouhaddou2023_gostres_filtered_KIN <- gost(
  query = query_filtered_KIN,
  organism = 'hsapiens',
  sources = c('KEGG')
)

baseline_KIN_results <- bouhaddou2023_gostres_baseline_KIN$result[which(bouhaddou2023_gostres_baseline_KIN$result$term_name %in% bouhaddou2023_gostres_filtered_KIN$result$term_name),]
filtered_KIN_results <- bouhaddou2023_gostres_filtered_KIN$result[which(bouhaddou2023_gostres_filtered_KIN$result$term_name %in% baseline_KIN_results$term_name),]

filtered_KIN_results <- filtered_KIN_results[1:5, ]
baseline_KIN_results <- baseline_KIN_results[which(baseline_KIN_results$term_name %in% filtered_KIN_results$term_name),]

filtered_KIN_results$color <- 'Filtered KIN'
baseline_KIN_results$color <- 'Baseline KIN'

sars_cov2_plot <- ggplot() +
  # First set of points
  geom_point(data = filtered_KIN_results, aes(x = reorder(str_wrap(term_name, width = 20), -p_value), y = -log10(p_value), color = color, size = intersection_size)) +
  
  # Second set of points
  geom_point(data = baseline_KIN_results, aes(x = reorder(str_wrap(term_name, width = 20), -p_value), y = -log10(p_value), color = color, size = intersection_size), shape = 17) +
  
  coord_flip() +
  labs(title = "SARS-CoV-2 (KEGG 2021)", 
       x = "Terms", 
       y = "-log10(p-value)",
       color = "Data set",
       size = str_wrap("Phosphorylation site count", width = 10)) +
  scale_color_manual(values = c("Filtered KIN" = colors[1], "Baseline KIN" = colors[2])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = "none")

# Wilkes 2015
wilkes2015_baseline_KIN <- read_tsv('../results/Hyperparameters/Wilkes2015/baseline_KIN/Wilkes2015_n15_alpha0.9_KIN.tsv')
wilkes2015_filtered_KIN <- read_tsv('../results/Hyperparameters/Wilkes2015/edge_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma1_beta0.4_PCST_filtered_KIN.tsv')

wilkes2015_gostres_baseline_KIN <- gost(
  query = unique(c(wilkes2015_baseline_KIN[which((wilkes2015_baseline_KIN$Edge_type == 'KS') & !(wilkes2015_baseline_KIN$Target_Uniprot %in% kinases)),]$Target_Uniprot)),
  organism = 'hsapiens',
  source = 'KEGG',
  significant = FALSE
)
wilkes2015_gostres_filtered_KIN <- gost(
  query = unique(c(wilkes2015_filtered_KIN[which((wilkes2015_filtered_KIN$Edge_type == 'KS') & !(wilkes2015_filtered_KIN$Target_Uniprot %in% kinases)),]$Target_Uniprot)),
  organism = 'hsapiens',
  source = 'KEGG',
  significant = FALSE
)

baseline_KIN_results <- wilkes2015_gostres_baseline_KIN$result[which(wilkes2015_gostres_baseline_KIN$result$term_name %in% wilkes2015_gostres_filtered_KIN$result$term_name),]
filtered_KIN_results <- wilkes2015_gostres_filtered_KIN$result[which(wilkes2015_gostres_filtered_KIN$result$term_name %in% baseline_KIN_results$term_name),]

filtered_KIN_results <- filtered_KIN_results[1:5, ]
baseline_KIN_results <- baseline_KIN_results[which(baseline_KIN_results$term_name %in% filtered_KIN_results$term_name),]

filtered_KIN_results$color <- 'Filtered KIN'
baseline_KIN_results$color <- 'Baseline KIN'

wilkes_plot <- ggplot() +
  # First set of points
  geom_point(data = filtered_KIN_results, aes(x = reorder(str_wrap(term_name, width = 20), -p_value), y = -log10(p_value), color = color, size = intersection_size)) +
  
  # Second set of points
  geom_point(data = baseline_KIN_results, aes(x = reorder(str_wrap(term_name, width = 20), -p_value), y = -log10(p_value), color = color, size = intersection_size), shape = 17) +
  
  coord_flip() +
  labs(title = "PI3K inhibition in resistant\n breast cancer cell line (KEGG 2021)", 
       x = "Terms", 
       y = "-log10(p-value)",
       color = "Data Set",
       size = str_wrap("Phosphorylation site count", width = 10)) +
  scale_color_manual(values = c("Filtered KIN" = colors[1], "Baseline KIN" = colors[2])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(color = "none")


### FINAL enrichment plot
combined_enrichment_plot <- grid.arrange(wilkes_plot, sars_cov2_plot, ncol=1)

# Showing kinase enrichment leads to significant results
serine_threonine_kinases <- read_tsv('../data/kinase_data/serine_threonine_kinases/kinase_name_mappings.tsv')
tyrosine_kinases <- read_tsv('../data/kinase_data/tyrosine_kinases/kinase_name_mappings.tsv')
gostres_kinases <- gost(
  query = c(serine_threonine_kinases$`ACC#`, tyrosine_kinases$`ACC#`),
  organism = 'hsapiens',
  source = 'KEGG'
)
kinase_enrichment_plot <- ggplot(head(gostres_kinases$result, 5), aes(x = reorder(term_name, -p_value), y = -log10(p_value), size = intersection_size)) +
  geom_point(color = "darkred") +
  coord_flip() +
  labs(title = "Kinase enrichments", 
       x = "Terms", 
       y = "-log10(p-value)", 
       size = "Gene Count") +
  theme_minimal()





