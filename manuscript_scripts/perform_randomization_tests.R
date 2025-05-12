library(tidyverse)
library(ComplexHeatmap)
library(igraph)
library(ggpubr)
library(gridExtra)
library(grid)
library(enrichR)
library(GGally)
source('./src/UniprotIDMapping.R')
# Randomization parameters
n_randomRandFiltering <- 1000

# Loading Omnipath network
op_net <- read_csv('./data/OmniPath/2023_10_11_KinaseDataOmniPath.csv', show_col_types = F)
op_net$Interaction <- paste0(op_net$enzyme, '_', op_net$TARGET_UP_ID, '_', op_net$TARGET_RES, op_net$TARGET_POS)

# Performs random filtering
perform_randomFiltering_test <- function(sub_scores, filtered_net, op_net, n_random = 1000) {
  sub_scores$Interaction <- paste0(sub_scores$Source, '_', sub_scores$Target)
  filtered_net$Interaction <- paste0(filtered_net$Source, '_', filtered_net$Target)
  filtered_net_interactions_in_op <- filtered_net[which(filtered_net$Interaction %in% op_net$Interaction),]
  no_randomInteractions <- c()
  for (i in 1:n_random) {
    net_random <- sub_scores[sample(nrow(sub_scores), nrow(filtered_net)), ]
    no_randomInteractions <- c(
      no_randomInteractions, 
      length(which(net_random$Interaction %in% op_net$Interaction))  
    ) 
  }
  
  mean_val <- mean(no_randomInteractions)
  sd_val <- sd(no_randomInteractions)
  z_score <- (nrow(filtered_net_interactions_in_op) - mean_val) / sd_val
  p_val <- pnorm(z_score, lower.tail = F)
  
  empirical_p_val <- p_val
  print(paste0('ECDF Val: ', ecdf(no_randomInteractions)(nrow(filtered_net_interactions_in_op))))
  
  return(list(no_randomInteractions = no_randomInteractions, p_val = p_val)) 
}


perform_randomRewiring_test <- function(filtered_net, op_net, n_random = 1000) {
  # preparing Dataframe
  # sub_scores$Pos <- sapply(strsplit(sub_scores$Target, '_'), function(x) { str_sub(x[[2]], 2) })
  # sub_scores$Target_w_Pos <- paste0(sub_scores$Target_LeadingProtein, '_', sub_scores$Pos)
  filtered_net$Interaction <- paste0(filtered_net$Source, '_', filtered_net$Target)
  filtered_net_interactions_in_op <- filtered_net[which(filtered_net$Interaction %in% op_net$Interaction),]
  
  g <- graph_from_data_frame(filtered_net[, c('Source', 'Target')])
  no_randomInteractions <- c()
  for (i in 1:n_random) {
    g_randomlyRewired <- g %>% rewire(keeping_degseq(niter = 10 * nrow(filtered_net))) 
    g_randomlyRewired_edges <- as.data.frame(as_edgelist(g_randomlyRewired))
    g_randomlyRewired_interactions <- paste0(g_randomlyRewired_edges$V1, '_', g_randomlyRewired_edges$V2) 
    no_randomInteractions <- c(
      no_randomInteractions, 
      length(which(g_randomlyRewired_interactions %in% op_net$Interaction))
    )
  }
  
  mean_val <- mean(no_randomInteractions)
  sd_val <- sd(no_randomInteractions)
  z_score <- (nrow(filtered_net_interactions_in_op) - mean_val) / sd_val
  p_val <- pnorm(z_score, lower.tail = F)
  
  empirical_p_val <- p_val
  print(paste0('ECDF Val: ', ecdf(no_randomInteractions)(nrow(filtered_net_interactions_in_op))))
  
  return(list(no_randomInteractions = no_randomInteractions, p_val = p_val)) 
}


compute_results <- function() {
  # get node filtered networks
  wilkes_baseline_KIN <- read_tsv('./results/Wilkes2015/baseline_KIN/MCF7G2_vs_Mock.tsv')
  bouhaddou_baseline_KIN <- read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/baseline_KIN/VIC10h_vs_Mock10h.tsv')

  wilkes_node_filtered_networks <- list(
    'DIFF' = read_tsv('./results/Wilkes2015/node_filtered_networks/MCF7G2_vs_Mock_DIFFnet.tsv'),
    'FS' = read_tsv('./results/Wilkes2015/node_filtered_networks/MCF7G2_vs_Mock_FSnet.tsv'),
    'FS + DIFF' = read_tsv('./results/Wilkes2015/node_filtered_networks/MCF7G2_vs_Mock_DIFFandFSnet.tsv')
  )

  bouhaddou_node_filtered_networks <- list(
    'DIFF' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/node_filtered_networks/VIC10h_vs_Mock10h_DIFFnet.tsv'),
    'FS' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/node_filtered_networks/VIC10h_vs_Mock10h_FSnet.tsv'),
    'FS + DIFF' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/node_filtered_networks/VIC10h_vs_Mock10h_DIFFandFSnet.tsv')
  )

  wilkes_node_filtered_networks <- lapply(wilkes_node_filtered_networks, function(x) {x[which(x$Type == 'KS'),]})
  bouhaddou_node_filtered_networks <- lapply(bouhaddou_node_filtered_networks, function(x) { x[which(x$Type == 'KS'),]})

  # perform randomization tests
  wilkes_randomFiltering_overlap <- c()
  wilkes_randomRewiring_overlap <- c()
  bouhaddou_randomFiltering_overlap <- c()
  bouhaddou_randomRewiring_overlap <- c()

  # perform randomization tests for Wilkes et al. 2015
  for (filtered_net in names(wilkes_node_filtered_networks)) {
    result_randomFiltering <- perform_randomFiltering_test(
      wilkes_baseline_KIN, 
      wilkes_node_filtered_networks[[filtered_net]], 
      op_net
    )
    wilkes_randomFiltering_overlap[[filtered_net]] <- result_randomFiltering[[1]]
    wilkes_randomFiltering_overlap[[paste0(filtered_net, '_pval')]] <- result_randomFiltering[[2]]
    
    result_randomRewiring <- perform_randomRewiring_test(
      wilkes_node_filtered_networks[[filtered_net]],
      op_net
    )
    wilkes_randomRewiring_overlap[[filtered_net]] <- result_randomRewiring[[1]]
    wilkes_randomRewiring_overlap[[paste0(filtered_net, '_pval')]] <- result_randomRewiring[[2]]
  }

  # perform randomization tests for Bouhaddou et al. 2023
  for (filtered_net in names(bouhaddou_node_filtered_networks)) {
    result_randomFiltering <- perform_randomFiltering_test(
      bouhaddou_baseline_KIN, 
      bouhaddou_node_filtered_networks[[filtered_net]], 
      op_net
    )
    bouhaddou_randomFiltering_overlap[[filtered_net]] <- result_randomFiltering[[1]]
    bouhaddou_randomFiltering_overlap[[paste0(filtered_net, '_pval')]] <- result_randomFiltering[[2]]
    
    result_randomRewiring <- perform_randomRewiring_test(
      bouhaddou_node_filtered_networks[[filtered_net]],
      op_net
    )
    bouhaddou_randomRewiring_overlap[[filtered_net]] <- result_randomRewiring[[1]]
    bouhaddou_randomRewiring_overlap[[paste0(filtered_net, '_pval')]] <- result_randomRewiring[[2]]
  }
}