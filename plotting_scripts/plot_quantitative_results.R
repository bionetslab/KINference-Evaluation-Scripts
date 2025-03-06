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


get_randDistributionPlot <- function(filtered_net, op_net, rand_test_overlaps, p_val, title = '', supp_plot = F) {
  filtered_net$Interaction <- paste0(filtered_net$Source, '_', filtered_net$Target)
  no_InteractionsInOP <- length(which(filtered_net$Interaction %in% op_net$Interaction))
  
  if (supp_plot) {
    y_text <- max(density(rand_test_overlaps)$y) / 2
    x_text_shift <- 0
    y_lim <- NA
  } else {
    y_lim <- ifelse(grepl('PI3K', title), 0.2, 0.1)
    y_text <- ifelse(grepl('PI3K', title), 0.15, 0.075)
    x_text_shift <- ifelse(grepl('PI3K', title), 7, 15)
  }
  
  plot <- ggplot(data.frame(value = rand_test_overlaps), aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), colour = 'black', fill = 'white') + 
    geom_density(lwd = 0.75, linetype = 1, alpha = 0.2, colour = '#E69F00') + 
    labs(x = 'Overlap', y = 'Occurence') +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_lim)) +
    geom_vline(aes(xintercept = no_InteractionsInOP), color = '#56B4E9') +
    geom_text(aes(x = no_InteractionsInOP - x_text_shift, y = y_text, label = no_InteractionsInOP), angle = -90, color = '#56B4E9', size = 8/.pt) +
    ggtitle(title) +
    theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    ) 
  
  return(plot)
}

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

## Publication plots: Plotting random background distribution plots
wilkes_plot_randFiltering <- get_randDistributionPlot(
  wilkes_node_filtered_networks[['FS + DIFF']],
  op_net,
  wilkes_randomFiltering_overlap[['FS + DIFF']],
  wilkes_randomFiltering_overlap[['FS + DIFF_pval']],
  title = 'PI3K inhibition in resistant\n breast cancer cell line\n FS + DIFF (RF)'
)

wilkes_plot_randRewiring <- get_randDistributionPlot(
  wilkes_node_filtered_networks[['FS + DIFF']],
  op_net,
  wilkes_randomRewiring_overlap[['FS + DIFF']],
  wilkes_randomRewiring_overlap[['FS + DIFF_pval']],
  title = 'PI3K inhibition in resistant\n breast cancer cell line\n FS + DIFF (RR)'
)

bouhaddou_plot_randFiltering <- get_randDistributionPlot(
  bouhaddou_node_filtered_networks[['FS + DIFF']],
  op_net,
  bouhaddou_randomFiltering_overlap[['FS + DIFF']],
  bouhaddou_randomFiltering_overlap[['FS + DIFF_pval']],
  title = 'SARS-CoV-2\n FS + DIFF (RF)'
)

bouhaddou_plot_randRewiring <- get_randDistributionPlot(
  bouhaddou_node_filtered_networks[['FS + DIFF']],
  op_net,
  bouhaddou_randomRewiring_overlap[['FS + DIFF']],
  bouhaddou_randomRewiring_overlap[['FS + DIFF_pval']],
  title = 'SARS-CoV-2\n FS + DIFF (RR)'
)
##

## Wilkes PCST
wilkes_PCST_networks <- list(
  'PCST' = read_tsv('./results/Wilkes2015/PCST_networks/MCF7G2_vs_Mock_PCST_net.tsv', show_col_types = F),
  'FS + PCST' = read_tsv('./results/Wilkes2015/PCST_networks/MCF7G2_vs_Mock_FSandPCST_net.tsv', show_col_types = F),
  'DIFF + PCST' = read_tsv('./results/Wilkes2015/PCST_networks/MCF7G2_vs_Mock_DIFFandPCST_net.tsv', show_col_types = F),
  'FS + DIFF + PCST' = read_tsv('./results/Wilkes2015/PCST_networks/MCF7G2_vs_Mock_DIFFandFSandPCST_net.tsv', show_col_types = F)
)
wilkes_PCST_networks <- lapply(wilkes_PCST_networks, function(x) { x[x$Type == 'KS',]})

wilkes_PCST_networks_randomFiltering_results <- list()
wilkes_PCST_networks_randomRewiring_results <- list()
for (filtered_net in names(wilkes_PCST_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    wilkes_baseline_KIN, 
    wilkes_PCST_networks[[filtered_net]], 
    op_net
  )
  wilkes_PCST_networks_randomFiltering_results[[filtered_net]] <- result_randomFiltering[[1]]
  wilkes_PCST_networks_randomFiltering_results[[paste0(filtered_net, '_pval')]] <- result_randomFiltering[[2]]
  
  result_randomRewiring <- perform_randomRewiring_test(
    wilkes_PCST_networks[[filtered_net]],
    op_net
  )
  wilkes_PCST_networks_randomRewiring_results[[filtered_net]] <- result_randomRewiring[[1]]
  wilkes_PCST_networks_randomRewiring_results[[paste0(filtered_net, '_pval')]] <- result_randomRewiring[[2]]
}

## Bouhaddou PCST
bouhaddou_PCST_networks <- list(
  'PCST' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/PCST_networks/VIC10h_vs_Mock10h_PCST_net.tsv', show_col_types = F),
  'FS + PCST' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/PCST_networks/VIC10h_vs_Mock10h_FSandPCST_net.tsv', show_col_types = F),
  'DIFF + PCST' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/PCST_networks/VIC10h_vs_Mock10h_DIFFandPCST_net.tsv', show_col_types = F),
  'FS + DIFF + PCST' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/PCST_networks/VIC10h_vs_Mock10h_DIFFandFSandPCST_net.tsv', show_col_types = F)
)

bouhaddou_PCST_networks <- lapply(bouhaddou_PCST_networks, function(x) { x[x$Type == 'KS',]})

bouhaddou_PCST_networks_randomFiltering_results <- list()
bouhaddou_PCST_networks_randomRewiring_results <- list()
for (filtered_net in names(bouhaddou_PCST_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    bouhaddou_baseline_KIN, 
    bouhaddou_PCST_networks[[filtered_net]], 
    op_net
  )
  bouhaddou_PCST_networks_randomFiltering_results[[filtered_net]] <- result_randomFiltering[[1]]
  bouhaddou_PCST_networks_randomFiltering_results[[paste0(filtered_net, '_pval')]] <- result_randomFiltering[[2]]
  
  result_randomRewiring <- perform_randomRewiring_test(
    bouhaddou_PCST_networks[[filtered_net]],
    op_net
  )
  bouhaddou_PCST_networks_randomRewiring_results[[filtered_net]] <- result_randomRewiring[[1]]
  bouhaddou_PCST_networks_randomRewiring_results[[paste0(filtered_net, '_pval')]] <- result_randomRewiring[[2]]
}

## bouhaddou CORR
bouhaddou_CORR_networks <- list(
  'CORR' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/correlation_networks/VIC10h_vs_Mock10h_X0.tsv'),
  'FS + CORR' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/correlation_networks/VIC10h_vs_Mock10h_X0_FSnet.tsv'),
  'DIFF + CORR' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/correlation_networks/VIC10h_vs_Mock10h_X0_DIFFnet.tsv'),
  'FS + DIFF + CORR' = read_tsv('./results/Bouhaddou2023/VIC10h_vs_Mock10h/correlation_networks/VIC10h_vs_Mock10h_X0_DIFFandFSnet.tsv')
)

bouhaddou_CORR_networks <- lapply(bouhaddou_CORR_networks, function(x) {x[x$Type == 'KS',]})

## Bouhaddou CORR
bouhaddou_corrFilteredNets_randomFiltering_results <- list()
bouhaddou_corrFilteredNets_randomRewiring_results <- list()
for (filtered_net in names(bouhaddou_CORR_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    bouhaddou_baseline_KIN, 
    bouhaddou_CORR_networks[[filtered_net]], 
    op_net
  )
  bouhaddou_corrFilteredNets_randomFiltering_results[[filtered_net]] <- result_randomFiltering[[1]]
  bouhaddou_corrFilteredNets_randomFiltering_results[[paste0(filtered_net, '_pval')]] <- result_randomFiltering[[2]]
  
  result_randomRewiring <- perform_randomRewiring_test(
    bouhaddou_CORR_networks[[filtered_net]],
    op_net
  )
  bouhaddou_corrFilteredNets_randomRewiring_results[[filtered_net]] <- result_randomRewiring[[1]]
  bouhaddou_corrFilteredNets_randomRewiring_results[[paste0(filtered_net, '_pval')]] <- result_randomRewiring[[2]]
}


## Publication plot: Heatmap
## Creating heatmap 
wilkes_edgeFiltering_randomFiltering_emp_p_vals <- c(
  wilkes_randomFiltering_overlap[['DIFF_pval']],
  wilkes_randomFiltering_overlap[['FS_pval']],
  wilkes_randomFiltering_overlap[['FS + DIFF_pval']],
  wilkes_PCST_networks_randomFiltering_results[['PCST_pval']],
  NA,
  #wilkes_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  #wilkes_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  wilkes_PCST_networks_randomFiltering_results[['FS + DIFF + PCST_pval']],
  
  #NA,
  #NA,
  NA
) 

wilkes_edgeFiltering_randomRewiring_emp_p_vals <- c(
  wilkes_randomRewiring_overlap[['DIFF_pval']],
  wilkes_randomRewiring_overlap[['FS_pval']],
  wilkes_randomRewiring_overlap[['FS + DIFF_pval']],
  wilkes_PCST_networks_randomRewiring_results[['PCST_pval']],
  NA,
  #wilkes_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  #wilkes_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  wilkes_PCST_networks_randomRewiring_results[['FS + DIFF + PCST_pval']],
  
  #NA,
  #NA,
  NA
)

bouhaddou_edgeFiltering_randomFiltering_emp_p_vals <- c(
  bouhaddou_randomFiltering_overlap[['DIFF_pval']],
  bouhaddou_randomFiltering_overlap[['FS_pval']],
  bouhaddou_randomFiltering_overlap[['FS + DIFF_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['PCST_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['CORR_pval']],
  #bouhaddou_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  #bouhaddou_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['FS + DIFF + PCST_pval']],
  #bouhaddou_corrFilteredNets_randomFiltering_results[['corrAndFS_net_pval']],
  #bouhaddou_corrFilteredNets_randomFiltering_results[['corrAndDiff_net_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['FS + DIFF + CORR_pval']]
)

bouhaddou_edgeFiltering_randomRewiring_emp_p_vals <- c(
  bouhaddou_randomRewiring_overlap[['DIFF_pval']],
  bouhaddou_randomRewiring_overlap[['FS_pval']],
  bouhaddou_randomRewiring_overlap[['FS + DIFF_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['PCST_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['CORR_pval']],
  #bouhaddou_PCST_networks_randomRewiring_results[['FS + PCST_pval']],
  #bouhaddou_PCST_networks_randomRewiring_results[['DIFF + PCST_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['FS + DIFF + PCST_pval']],
  #bouhaddou_corrFilteredNets_randomRewiring_results[['corrAndFS_net_pval']],
  #bouhaddou_corrFilteredNets_randomRewiring_results[['corrAndDiff_net_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['FS + DIFF + CORR_pval']]
)

edgeFiltering_emp_p_vals <- data.frame(
  'PI3K inhibition in resistant\n breast cancer cell line' = wilkes_edgeFiltering_randomFiltering_emp_p_vals,
  'SARS-CoV-2' = bouhaddou_edgeFiltering_randomFiltering_emp_p_vals,
  'PI3K inhibition in resistant\n breast cancer cell line' = wilkes_edgeFiltering_randomRewiring_emp_p_vals,
  'SARS-CoV-2' = bouhaddou_edgeFiltering_randomRewiring_emp_p_vals,
  check.names = F
)

rownames(edgeFiltering_emp_p_vals) <- c('DIFF', 'FS', 'FS + DIFF', 'PCST', 'CORR', 'FS + DIFF + PCST',  'FS + DIFF + CORR')

cell_annotation <- data.frame(
  rf1 = c('0.078', '< 0.001', '< 0.001', '0.435', 'NA', '< 0.001', 'NA'), 
  rf2 = c('0.048', '< 0.001', '< 0.001', '0.584', '< 0.001', '< 0.001', '< 0.001'),
  rr1 = c('< 0.001', '< 0.001', '< 0.001', '< 0.001', 'NA', '< 0.001', 'NA'),
  rr2 = c('< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001','< 0.001', '< 0.001')
)

hm_data <- t(edgeFiltering_emp_p_vals)

hm <- ComplexHeatmap::pheatmap(
  as.matrix(hm_data), 
  color = colorRampPalette(c("navy", "", "red"))(500),
  breaks = seq(0, 0.1, length.out = 501),
  cellwidth = 30,
  cellheight = 30,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    text_color <- ifelse(is.na(edgeFiltering_emp_p_vals[j, i]), "black", "white")
    grid.text(cell_annotation[j, i], x, y, gp = gpar(col = text_color, fontsize = 8))
  },
  heatmap_legend_param = list(
    title = 'P-value', at = c(0, 0.025, 0.05, 0.075, 0.1), labels = c('0', '0.025', '0.05', '0.075', '> 0.1')
  ),
  row_split = c('RF', 'RF', 'RR', 'RR'),
  column_split = factor(c(rep('Node filters', 3), rep('Edge filters\n without node\n filtering', 2), rep('Edge filters\n with node\n filtering', 2)), levels=unique(c(rep('Node filters', 3), rep('Edge filters\n without node\n filtering', 2), rep('Edge filters\n with node\n filtering', 2)))),
  fontsize = 7,
  row_title_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 7),
  angle_col = c('45')
)

## Publication table:
compute_reduction <- function(x, y) {
  return (paste0(round(100 * (1 - y / x), 2), '%'))
}

wilkes_nrows <- lapply(wilkes_node_filtered_networks, function(x) { nrow(x) })
bouhaddou_nrows <- lapply(bouhaddou_node_filtered_networks, function(x) { nrow(x) })

bouhaddou_CORR_nrows <- lapply(bouhaddou_CORR_networks, function(x) { nrow(x) })
wilkes_PCST_nrows <- lapply(wilkes_PCST_networks, function(x) { nrow(x) })
bouhaddou_PCST_nrows <- lapply(bouhaddou_PCST_networks, function(x) { nrow(x) })


wilkes_baseline_rows <- nrow(wilkes_baseline_KIN)
bouhaddou_baseline_rows <- nrow(bouhaddou_baseline_KIN)

data <- data.frame(matrix(c(
  wilkes_baseline_rows, "0%", bouhaddou_baseline_rows, "0%",
  wilkes_nrows[['FS']],  compute_reduction(wilkes_baseline_rows, wilkes_nrows[['FS']]), bouhaddou_nrows[['FS']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_nrows[['FS']]),
  wilkes_nrows[['DIFF']],  compute_reduction(wilkes_baseline_rows, wilkes_nrows[['DIFF']]), bouhaddou_nrows[['DIFF']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_nrows[['DIFF']]),
  wilkes_nrows[['FS + DIFF']], compute_reduction(wilkes_baseline_rows, wilkes_nrows[['FS + DIFF']]), bouhaddou_nrows[['FS + DIFF']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_nrows[['FS + DIFF']]),
  "NA", "NA", bouhaddou_CORR_nrows[['CORR']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_CORR_nrows[['CORR']]),
  "NA", "NA", bouhaddou_CORR_nrows[['FS + CORR']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_CORR_nrows[['FS + CORR']]),
  "NA", "NA", bouhaddou_CORR_nrows[['DIFF + CORR']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_CORR_nrows[['DIFF + CORR']]),
  "NA", "NA", bouhaddou_CORR_nrows[['FS + DIFF + CORR']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_CORR_nrows[['FS + DIFF + CORR']]),
  wilkes_PCST_nrows[['PCST']], compute_reduction(wilkes_baseline_rows, wilkes_PCST_nrows[['PCST']]), bouhaddou_PCST_nrows[['PCST']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_PCST_nrows[['PCST']]),
  wilkes_PCST_nrows[['FS + PCST']], compute_reduction(wilkes_baseline_rows, wilkes_PCST_nrows[['FS + PCST']]), bouhaddou_PCST_nrows[['FS + PCST']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_PCST_nrows[['FS + PCST']]),
  wilkes_PCST_nrows[['DIFF + PCST']], compute_reduction(wilkes_baseline_rows, wilkes_PCST_nrows[['DIFF + PCST']]), bouhaddou_PCST_nrows[['DIFF + PCST']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_PCST_nrows[['DIFF + PCST']]),
  wilkes_PCST_nrows[['FS + DIFF + PCST']], compute_reduction(wilkes_baseline_rows, wilkes_PCST_nrows[['FS + DIFF + PCST']]), bouhaddou_PCST_nrows[['FS + DIFF + PCST']], compute_reduction(bouhaddou_baseline_rows, bouhaddou_PCST_nrows[['FS + DIFF + PCST']])
), ncol = 4, byrow = TRUE))

rownames(data) <- c('Baseline KIN', 'FS', 'DIFF', 'FS + DIFF', 'CORR', 'FS + CORR', 'DIFF + CORR', 'FS + DIFF + CORR', 'PCST', 'FS + PCST', 'DIFF + PCST', 'FS + DIFF + PCST')
colnames(data) <- c('PI3K inhibition in resistant\n breast cancer cell line', 'Reduction \n(%)', 'SARS-CoV-2', 'Reduction \n(%)')

table <- tableGrob(
  data,
  theme = ttheme_default(base_size = 8)
)

ggarrange(table, ncol = 1, nrow = 1)


# Arrange the plots
results_plot <- ggarrange(
  ggarrange(
    grid.grabExpr(draw(hm)),
    ncol = 1, nrow = 1,
    labels = c('A')
  ),
  ggarrange(
    wilkes_plot_randFiltering, wilkes_plot_randRewiring, bouhaddou_plot_randFiltering, bouhaddou_plot_randRewiring,
    ncol = 2, nrow = 2, 
    labels = c('B', '', 'C', '')
  ),
  ncol = 2, nrow = 1, widths = c(1.4, 1)
)

##

## GSEA Analysis (Publication plot)
bouhaddou_PCST_FSandDiffFiltered_net_genes <- translateUniprot2GeneName(
  unique(c(bouhaddou_PCST_networks[['FS + DIFF + PCST']]$Source, bouhaddou_PCST_networks[['FS + DIFF + PCST']]$Target_Uniprot)),
  species = 'HUMAN'
)
write_tsv(data.frame(genes = bouhaddou_PCST_networks[['FS + DIFF + PCST']]), './bouhaddou_PCST_FSandDiffFiltered_net_genes_onlyTargets.tsv')

wilkes_PCST_FSandDiffFiltered_net_genes <- translateUniprot2GeneName(
  unique(c(wilkes_PCST_networks[['FS + DIFF + PCST']]$Source, wilkes_PCST_networks[['FS + DIFF + PCST']]$Target_Uniprot)),
  species = 'HUMAN'
)
write_tsv(data.frame(genes = wilkes_PCST_networks[['FS + DIFF + PCST']]), './wilkes_PCST_FSandDiffFiltered_net_genes.tsv')

comp_bouhaddou_genes <- translateUniprot2GeneName(
  unique(c(bouhaddou_baseline_KIN$Source, bouhaddou_baseline_KIN$Target_Uniprot)),
  species = 'HUMAN'
)
write_tsv(data.frame(genes = comp_bouhaddou_genes), './bouhaddou_subScores_genes.tsv')

comp_wilkes_genes <- translateUniprot2GeneName(
  unique(c(wilkes_baseline_KIN$Source, wilkes_baseline_KIN$Target_Uniprot)),
  species = 'HUMAN'
)
write_tsv(data.frame(genes = comp_wilkes_genes), './wilkes_subScores_genes.tsv')

# perform enrichment analysis
websiteLive <- getOption('enrichR.live')
if (websiteLive) dbs <- listEnrichrDbs()
dbs <- c("KEGG_2021_Human")
if (websiteLive) {
  bouhaddou_enriched_PCSTandFSandDIFF <- enrichr(bouhaddou_PCST_FSandDiffFiltered_net_genes, dbs)
  wilkes_enriched_PCSTandFSandDIFF <- enrichr(wilkes_PCST_FSandDiffFiltered_net_genes, dbs)
  bouhaddou_enriched_comp <- enrichr(comp_bouhaddou_genes, dbs)
  wilkes_enriched_comp <- enrichr(comp_wilkes_genes, dbs)
}

# Plotting
n_terms <- 5

wilkes_enriched_PCSTandFSandDIFF[[1]] <- wilkes_enriched_PCSTandFSandDIFF[[1]][1:n_terms,]
bouhaddou_enriched_PCSTandFSandDIFF[[1]] <- bouhaddou_enriched_PCSTandFSandDIFF[[1]][1:n_terms,]

wilkes_enriched_PCSTandFSandDIFF[[1]]$Count <- sapply(wilkes_enriched_PCSTandFSandDIFF[[1]]$Genes, function(x) {str_count(x, ';') + 1})
bouhaddou_enriched_PCSTandFSandDIFF[[1]]$Count <- sapply(bouhaddou_enriched_PCSTandFSandDIFF[[1]]$Genes, function(x) {str_count(x, ';') + 1})

wilkes_enriched_PCSTandFSandDIFF[[1]]$Term <- factor(wilkes_enriched_PCSTandFSandDIFF[[1]]$Term, levels = wilkes_enriched_PCSTandFSandDIFF[[1]]$Term)
bouhaddou_enriched_PCSTandFSandDIFF[[1]]$Term <- factor(bouhaddou_enriched_PCSTandFSandDIFF[[1]]$Term, levels = bouhaddou_enriched_PCSTandFSandDIFF[[1]]$Term)

enrich_plot1 <- ggplot(data = wilkes_enriched_PCSTandFSandDIFF[[1]][1:n_terms,], aes(x=reorder(Term, P.value), y=Count, fill = P.value)) +
  geom_bar(stat='identity', width = 0.8) +
  scale_fill_gradient(low = "red", high = "blue") +
  guides(fill = guide_colorbar(reverse = TRUE)) +
  scale_x_discrete(
    limits = rev(levels(wilkes_enriched_PCSTandFSandDIFF[[1]]$Term)),
    labels = rev(c('Neurotrophin signaling pathway', 'Progesterone-mediated\n oocyte maturation', 'T cell receptor signaling pathway', 'Fc epsilon RI\n signaling pathway', 'Prolactin signaling pathway'))
    ) + 
  theme_minimal() +
  labs(
    title = 'PI3K inhibition in resistant\n breast cancer cell line (KEGG 2021)',
    x = '',
    y = 'Gene count',
    fill = 'P-value'
  ) +
  theme(
    legend.text = element_text(size=6),
    axis.text.x = element_text(size=6),
    legend.key.size = unit(0.25, 'cm'),
    axis.text.y = element_text(size=6),
    legend.title = element_text(size=6),
    axis.title.x = element_text(size=6),
    plot.title = element_text(size=6, hjust = 0.5)
  ) +
  coord_flip()


enrich_plot2 <- ggplot(data = bouhaddou_enriched_PCSTandFSandDIFF[[1]][1:n_terms,], aes(x=reorder(Term, P.value), y=Count, fill = P.value)) +
  geom_bar(stat='identity', width = 0.8) +
  scale_fill_gradient(low = "red", high = "blue") +
  guides(fill = guide_colorbar(reverse = TRUE)) +
  scale_x_discrete(
    limits = rev(levels(bouhaddou_enriched_PCSTandFSandDIFF[[1]]$Term)),
    labels = c('MAPK signaling pathway', 'Autophagy', 'Neurotrophin signaling pathway', 'FoxO signaling pathway', 'ErbB signaling pathway')
    ) + 
  theme_minimal() +
  labs(
    title = 'SARS-CoV-2 (KEGG 2021)',
    x = '',
    y = 'Gene count',
    fill = 'P-value'
  ) +
  theme(
    legend.text = element_text(size=7),
    legend.key.size = unit(0.25, 'cm'),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7),
    legend.title = element_text(size=7),
    axis.title.x = element_text(size=7),
    plot.title = element_text(size=7, hjust = 0.5)
  ) +
  coord_flip()

p_enrich <- ggarrange(enrich_plot1, enrich_plot2, 
                      #labels = c('A', 'B'),
                      ncol = 1, nrow = 2, font.label=list(color="black",size=8))

## Supplementary plots: all background distributions
wilkes_RF_plots <- list()
wilkes_RR_plots <- list()
for (name in names(wilkes_node_filtered_networks)) {
  wilkes_RF_plots[[name]] <- get_randDistributionPlot(
    wilkes_node_filtered_networks[[name]],
    op_net,
    wilkes_randomFiltering_overlap[[name]],
    wilkes_randomFiltering_overlap[[paste0(name, '_pval')]],
    title = paste0('PI3K inhibition in resistant\n breast cancer cell, ', str_to_upper(str_replace_all(str_replace_all(name, '_net', ''), 'And', ' + ')), ' (RF)'),
    supp_plot = T
  )
  wilkes_RR_plots[[name]] <- get_randDistributionPlot(
    wilkes_node_filtered_networks[[name]],
    op_net,
    wilkes_randomRewiring_overlap[[name]],
    wilkes_randomRewiring_overlap[[paste0(name, '_pval')]],
    title = paste0('PI3K inhibition in resistant\n breast cancer cell, ', str_to_upper(str_replace_all(str_replace_all(name, '_net', ''), 'And', ' + ')), ' (RR)'),
    supp_plot = T
  )
}

bouhaddou_RF_plots <- list()
bouhaddou_RR_plots <- list()
# pcst nets
for (name in names(bouhaddou_node_filtered_networks)) {
  bouhaddou_RF_plots[[name]] <- get_randDistributionPlot(
    bouhaddou_node_filtered_networks[[name]],
    op_net,
    bouhaddou_randomFiltering_overlap[[name]],
    bouhaddou_randomFiltering_overlap[[paste0(name, '_pval')]],
    title = paste0('SARS-CoV-2, ', str_to_upper(str_replace_all(str_replace_all(name, '_net', ''), 'And', ' + ')), ' (RF)'),
    supp_plot = T
  )
  bouhaddou_RR_plots[[name]] <- get_randDistributionPlot(
    bouhaddou_node_filtered_networks[[name]],
    op_net,
    bouhaddou_randomRewiring_overlap[[name]],
    bouhaddou_randomRewiring_overlap[[paste0(name, '_pval')]],
    title = paste0('SARS-CoV-2, ', str_to_upper(str_replace_all(str_replace_all(name, '_net', ''), 'And', ' + ')), ' (RR)'),
    supp_plot = T
  )
}

# Supp Figure 01
bouhaddou_nodeFiltering_results_plot <- ggarrange(
  ggarrange(
    bouhaddou_RF_plots[['FS']], bouhaddou_RR_plots[['FS']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_RF_plots[['DIFF']], bouhaddou_RR_plots[['DIFF']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_RF_plots[['FS + DIFF']], bouhaddou_RR_plots[['FS + DIFF']],
    ncol = 2, nrow = 1
  ),
  ncol = 1, nrow = 3
)

# Supp Figure 02
wilkes_nodeFiltering_results_plot <- ggarrange(
  ggarrange(
    wilkes_RF_plots[['FS']], wilkes_RR_plots[['FS']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    wilkes_RF_plots[['DIFF']], wilkes_RR_plots[['DIFF']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    wilkes_RF_plots[['FS + DIFF']], wilkes_RR_plots[['FS + DIFF']],
    ncol = 2, nrow = 1
  ),
  ncol = 1, nrow = 3
)
##

## Supplementary Plot: Edge filtering background distribution plots
wilkes_PCST_RF_plots <- list()
wilkes_PCST_RR_plots <- list()

for (name in names(wilkes_PCST_networks)) {
  wilkes_PCST_RF_plots[[name]] <- get_randDistributionPlot(
    wilkes_PCST_networks[[name]],
    op_net,
    wilkes_PCST_networks_randomFiltering_results[[name]],
    wilkes_PCST_networks_randomFiltering_results[[paste0(name, '_pval')]],
    title = paste0('PI3K inhibition in resistant\n breast cancer cell, ', name, ' (RF)'),
    supp_plot = T
  )
  wilkes_PCST_RR_plots[[name]] <- get_randDistributionPlot(
    wilkes_PCST_networks[[name]],
    op_net,
    wilkes_PCST_networks_randomRewiring_results[[name]],
    wilkes_PCST_networks_randomRewiring_results[[paste0(name, '_pval')]],
    title = paste0('PI3K inhibition in resistant\n breast cancer cell, ', name, ' (RR)'),
    supp_plot = T
  )
}

bouhaddou_EdgeFilter_RF_plots <- list()
bouhaddou_EdgeFilter_RR_plots <- list()
# pcst nets
for (name in names(bouhaddou_PCST_networks)) {
  bouhaddou_EdgeFilter_RF_plots[[name]] <- get_randDistributionPlot(
    bouhaddou_PCST_networks[[name]],
    op_net,
    bouhaddou_PCST_networks_randomFiltering_results[[name]],
    bouhaddou_PCST_networks_randomFiltering_results[[paste0(name, '_pval')]],
    title = paste0('SARS-CoV-2, ', name, ' (RF)'),
    supp_plot = T
  )
  bouhaddou_EdgeFilter_RR_plots[[name]] <- get_randDistributionPlot(
    bouhaddou_PCST_networks[[name]],
    op_net,
    bouhaddou_PCST_networks_randomRewiring_results[[name]],
    bouhaddou_PCST_networks_randomRewiring_results[[paste0(name, '_pval')]],
    title = paste0('SARS-CoV-2, ', name, ' (RR)'),
    supp_plot = T
  )
}
# corr nets
for (name in names(bouhaddou_CORR_networks)) {
  bouhaddou_EdgeFilter_RF_plots[[name]] <- get_randDistributionPlot(
    bouhaddou_CORR_networks[[name]],
    op_net,
    bouhaddou_corrFilteredNets_randomFiltering_results[[name]],
    bouhaddou_corrFilteredNets_randomFiltering_results[[paste0(name, '_pval')]],
    title = paste0('SARS-CoV-2, ', str_to_upper(str_replace_all(str_replace_all(name, '_net', ''), 'And', ' + ')), ' (RF)'),
    supp_plot = T
  )
  bouhaddou_EdgeFilter_RR_plots[[name]] <- get_randDistributionPlot(
    bouhaddou_CORR_networks[[name]],
    op_net,
    bouhaddou_corrFilteredNets_randomRewiring_results[[name]],
    bouhaddou_corrFilteredNets_randomRewiring_results[[paste0(name, '_pval')]],
    title = paste0('SARS-CoV-2, ', str_to_upper(str_replace_all(str_replace_all(name, '_net', ''), 'And', ' + ')), ' (RR)'),
    supp_plot = T
  )
}

# supp. Fig. 05
wilkes_PCST_results_plot <- ggarrange(
  ggarrange(
    wilkes_PCST_RF_plots[['PCST']], wilkes_PCST_RR_plots[['PCST']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    wilkes_PCST_RF_plots[['DIFF + PCST']], wilkes_PCST_RR_plots[['DIFF + PCST']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    wilkes_PCST_RF_plots[['FS + PCST']], wilkes_PCST_RR_plots[['FS + PCST']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    wilkes_PCST_RF_plots[['FS + DIFF + PCST']], wilkes_PCST_RR_plots[['FS + DIFF + PCST']],
    ncol = 2, nrow = 1
  ),
  ncol = 1, nrow = 4
)

# supp. Fig. 03
bouhaddou_PCST_results_plot <- ggarrange(
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['PCST']], bouhaddou_EdgeFilter_RR_plots[['PCST']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['DIFF + PCST']], bouhaddou_EdgeFilter_RR_plots[['DIFF + PCST']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['FS + PCST']], bouhaddou_EdgeFilter_RR_plots[['FS + PCST']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['FS + DIFF + PCST']], bouhaddou_EdgeFilter_RR_plots[['FS + DIFF + PCST']],
    ncol = 2, nrow = 1
  ),
  ncol = 1, nrow = 4
)

# supp. Fig. 04
bouhaddou_CORR_results_plot <- ggarrange(
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['CORR']], bouhaddou_EdgeFilter_RR_plots[['CORR']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['FS + CORR']], bouhaddou_EdgeFilter_RR_plots[['FS + CORR']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['DIFF + CORR']], bouhaddou_EdgeFilter_RR_plots[['DIFF + CORR']],
    ncol = 2, nrow = 1
  ),
  ggarrange(
    bouhaddou_EdgeFilter_RF_plots[['FS + DIFF + CORR']], bouhaddou_EdgeFilter_RR_plots[['FS + DIFF + CORR']],
    ncol = 2, nrow = 1
  ),
  ncol = 1, nrow = 4
)

## Creating heatmap 
wilkes_edgeFiltering_randomFiltering_emp_p_vals <- c(
  wilkes_randomFiltering_overlap[['DIFF_pval']],
  wilkes_randomFiltering_overlap[['FS_pval']],
  wilkes_randomFiltering_overlap[['FS + DIFF_pval']],
  wilkes_PCST_networks_randomFiltering_results[['PCST_pval']],
  wilkes_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  wilkes_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  wilkes_PCST_networks_randomFiltering_results[['FS + DIFF + PCST_pval']],
  NA,
  NA,
  NA,
  NA
) 

wilkes_edgeFiltering_randomRewiring_emp_p_vals <- c(
  wilkes_randomRewiring_overlap[['DIFF_pval']],
  wilkes_randomRewiring_overlap[['FS_pval']],
  wilkes_randomRewiring_overlap[['FS + DIFF_pval']],
  wilkes_PCST_networks_randomRewiring_results[['PCST_pval']],
  wilkes_PCST_networks_randomRewiring_results[['FS + PCST_pval']],
  wilkes_PCST_networks_randomRewiring_results[['DIFF + PCST_pval']],
  wilkes_PCST_networks_randomRewiring_results[['FS + DIFF + PCST_pval']],
  NA,
  NA,
  NA,
  NA
)

bouhaddou_edgeFiltering_randomFiltering_emp_p_vals <- c(
  bouhaddou_randomFiltering_overlap[['DIFF_pval']],
  bouhaddou_randomFiltering_overlap[['FS_pval']],
  bouhaddou_randomFiltering_overlap[['FS + DIFF_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['PCST_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['FS + DIFF + PCST_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['CORR_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['FS + CORR_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['DIFF + CORR_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['FS + DIFF + CORR_pval']]
)

bouhaddou_edgeFiltering_randomRewiring_emp_p_vals <- c(
  bouhaddou_randomRewiring_overlap[['DIFF_pval']],
  bouhaddou_randomRewiring_overlap[['FS_pval']],
  bouhaddou_randomRewiring_overlap[['FS + DIFF_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['PCST_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['FS + PCST_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['DIFF + PCST_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['FS + DIFF + PCST_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['CORR_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['FS + CORR_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['DIFF + CORR_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['FS + DIFF + CORR_pval']]
)

edgeFiltering_emp_p_vals <- data.frame(
  'PI3K inhibition in resistant\n breast cancer cell' = wilkes_edgeFiltering_randomFiltering_emp_p_vals,
  'SARS-CoV-2' = bouhaddou_edgeFiltering_randomFiltering_emp_p_vals,
  'PI3K inhibition in resistant\n breast cancer cell' = wilkes_edgeFiltering_randomRewiring_emp_p_vals,
  'SARS-CoV-2' = bouhaddou_edgeFiltering_randomRewiring_emp_p_vals,
  check.names = F
)

rownames(edgeFiltering_emp_p_vals) <- c('DIFF', 'FS', 'FS + DIFF', 'PCST', 'FS + PCST', 'DIFF + PCST', 'FS + DIFF + PCST', 'CORR', 'FS + CORR', 'DIFF + CORR', 'FS + DIFF + CORR')

cell_annotation <- data.frame(
  rf1 = c('0.076', '< 0.001', '< 0.001', '0.436', '< 0.001', '0.359', '< 0.001', 'NA', 'NA', 'NA', 'NA'), 
  rf2 = c('0.040', '< 0.001', '< 0.001', '0.577', '< 0.001', '0.121', '< 0.001', '< 0.001', '< 0.001', '0.999', '< 0.001'),
  rr1 = c('< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', 'NA', 'NA', 'NA', 'NA'),
  rr2 = c('< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001', '< 0.001')
)

# Supp. Fig. 06
hm <- ComplexHeatmap::pheatmap(
  edgeFiltering_emp_p_vals, 
  color = colorRampPalette(c("navy", "white", "red"))(500),
  breaks = seq(0, 0.1, length.out = 501),
  cellwidth = 30,
  cellheight = 30,
  cluster_rows = F,
  cluster_cols = F,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    text_color <- ifelse(is.na(edgeFiltering_emp_p_vals[i, j]), "black", "white")
    grid.text(cell_annotation[i, j], x, y, gp = gpar(col = text_color, fontsize = 8))
  },
  heatmap_legend_param = list(
    title = 'P-value', at = c(0, 0.025, 0.05, 0.075, 0.1)
  ),
  column_split = c('RF', 'RF', 'RR', 'RR'),
  row_split = factor(c(rep('Node filters', 3), rep('Edge filters (PCST)', 4), rep('Edge filters (CORR)', 4)), levels=unique(c(rep('Node filters', 3), rep('Edge filters (PCST)', 4), rep('Edge filters (CORR)', 4)))),
  fontsize = 8,
  row_title_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 8),
)


### graph generation
csnk2a1_sub_net <- bouhaddou_PCST_networks$`FS + DIFF + PCST`[which(bouhaddou_PCST_networks$`FS + DIFF + PCST`$Source == 'P68400'),]
csnk2a1_sub_net <- csnk2a1_sub_net[order(csnk2a1_sub_net$f, decreasing = T),]
# top 10 phosphorylated edges
csnk2a1_sub_net <- csnk2a1_sub_net[1:10, ]
csnk2a1_sub_net$Source <- translateUniprot2GeneName(csnk2a1_sub_net$Source)
csnk2a1_sub_net$Target_Uniprot <- translateUniprot2GeneName(csnk2a1_sub_net$Target_Uniprot)
csnk2a1_sub_net$Pos <- sapply(strsplit(csnk2a1_sub_net$Target, '_'), function(x) { x[2] })

g_data <- csnk2a1_sub_net[, c('Source', 'Target_Uniprot')]

col_palette <-RColorBrewer::brewer.pal(name = 'Set2', n = 4)
g <- graph_from_data_frame(g_data, directed = TRUE)
ggnet2(g, size = c(5, csnk2a1_sub_net$f), label = TRUE, color = c(col_palette[1], rep(col_palette[4], 10)), shape = c(17, rep(19, 10)), edge.label = csnk2a1_sub_net$Pos, legend.position = 'none', arrow.size=8, arrow.gap=0.0525)

