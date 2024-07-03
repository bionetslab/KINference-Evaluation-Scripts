library(tidyverse)
library(ComplexHeatmap)
library(igraph)
library(ggpubr)
# Randomization parameters
n_randomRandFiltering <- 1000

# Loading Omnipath network
op_net <- read_csv('./data/OmniPath/2023_10_11_KinaseDataOmniPath.csv', show_col_types = F)
op_net$Interaction <- paste0(op_net$enzyme, '_', op_net$TARGET_UP_ID, '_', op_net$TARGET_RES, op_net$TARGET_POS)

# get all node filtered resulting networks
retrieve_NodeFiltered_networks <- function(sub_scores, functional_scores, beta = NA, gamma = NA) {
  
  diff_net <- NULL
  fs_net <- NULL
  fsAndDiff_net <- NULL
  
  if (is.na(beta) && is.na(gamma)) {
    stop("At least beta or gamma have to be given to retrieve a filtered network!")
  } else if(is.na(beta)) {
    diff_net <- sub_scores[abs(sub_scores$log2FC) >= gamma, ]
  } else {
    sub_scores$Pos <- sapply(
      strsplit(sub_scores$Target, '_'), 
      function(x) {
        str_sub(x[[2]], 2)    
      }
    )
    sub_scores$Target_wPos <- paste0(sub_scores$Target_LeadingProtein, '_', sub_scores$Pos)
    
    functional_scores <- column_to_rownames(functional_scores, var = 'sites')
    sub_scores$FS <- functional_scores[sub_scores$Target_wPos,]
    
    fs_net <- sub_scores[sub_scores$FS >= beta, ]
    if (!is.na(gamma)) {
      diff_net <- sub_scores[abs(sub_scores$log2FC) >= gamma, ]
      fsAndDiff_net <- fs_net[fs_net$log2FC >= gamma, ]
    }
  }
  
  return (list(diff_net = diff_net, fs_net = fs_net, fsAndDiff_net = fsAndDiff_net))
}

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
      length(which(g_randomlyRewired_interactions %in% op_network$Interaction))
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


get_randDistributionPlot <- function(filtered_net, op_net, rand_test_overlaps, p_val, title = '') {
  filtered_net$Interaction <- paste0(filtered_net$Source, '_', filtered_net$Target)
  no_InteractionsInOP <- length(which(filtered_net$Interaction %in% op_net$Interaction))
  
  y_lim <- ifelse(grepl('Wilkes', title), 0.2, 0.1)
  y_text <- ifelse(grepl('Wilke', title), 0.15, 0.075)
  x_text_shift <- ifelse(grepl('Wilke', title), 1, 6)
  plot <- ggplot(data.frame(value = rand_test_overlaps), aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), colour = 'black', fill = 'white') + 
    geom_density(lwd = 0.75, linetype = 1, alpha = 0.2, colour = '#E69F00') + 
    labs(x = 'Overlap', y = 'Occurence') +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, y_lim)) +
    geom_vline(aes(xintercept = no_InteractionsInOP), color = '#56B4E9') +
    geom_text(aes(x = no_InteractionsInOP - x_text_shift, y = y_text, label = no_InteractionsInOP), angle = -90, color = '#56B4E9') +
    ggtitle(title) +
    theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    ) 
    
  
  return(plot)
}


### Plotting
# Path and parameter definitions
wilkes_sub_scores <- read_tsv('./results/Wilkes2015/substrate_scores/MCF7-G2_vs_Mock.tsv', show_col_types = F)
wilkes_functional_scores <- read_tsv('./results/Wilkes2015/functionalScores.tsv', show_col_types = F)
wilkes_beta <- 0.4
wilkes_gamma <- 1.0

bouhaddou_sub_scores <- read_tsv('./results/Bouhaddou2023/substrate_scores/VIC_10h_vs_Mock_10h.tsv', show_col_types = F)
bouhaddou_functional_scores <- read_tsv('./results/Bouhaddou2023/functionalScores.tsv', show_col_types = F)
bouhaddou_beta <- 0.4
bouhaddou_gamma <- 1.0

# get node filtered networks
wilkes_node_filtered_networks <- retrieve_NodeFiltered_networks(
  wilkes_sub_scores,
  wilkes_functional_scores,
  wilkes_beta,
  wilkes_gamma
)

bouhaddou_node_filtered_networks <- retrieve_NodeFiltered_networks(
  bouhaddou_sub_scores,
  bouhaddou_functional_scores,
  bouhaddou_beta,
  bouhaddou_gamma
)

# perform randomization tests
wilkes_randomFiltering_overlap <- c()
wilkes_randomRewiring_overlap <- c()
bouhaddou_randomFiltering_overlap <- c()
bouhaddou_randomRewiring_overlap <- c()

# perform randomization tests for Wilkes et al. 2015
for (filtered_net in names(wilkes_node_filtered_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    wilkes_sub_scores, 
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
    bouhaddou_sub_scores, 
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

## Plotting random distribution plots

wilkes_plot_randFiltering <- get_randDistributionPlot(
  wilkes_node_filtered_networks$fsAndDiff_net,
  op_net,
  wilkes_randomFiltering_overlap$fsAndDiff_net,
  wilkes_randomFiltering_overlap$fsAndDiff_net_pval,
  title = 'Wilkes et al. (Random Filtering)'
)

wilkes_plot_randRewiring <- get_randDistributionPlot(
  wilkes_node_filtered_networks$fsAndDiff_net,
  op_net,
  wilkes_randomRewiring_overlap$fsAndDiff_net,
  wilkes_randomRewiring_overlap$fsAndDiff_net_pval,
  title = 'Wilkes et al. (Random Rewiring)'
)

bouhaddou_plot_randFiltering <- get_randDistributionPlot(
  bouhaddou_node_filtered_networks$fsAndDiff_net,
  op_net,
  bouhaddou_randomFiltering_overlap$fsAndDiff_net,
  bouhaddou_randomFiltering_overlap$fsAndDiff_net_pval,
  title = 'Bouhaddou et al. (Random Filtering)'
)

bouhaddou_plot_randRewiring <- get_randDistributionPlot(
  bouhaddou_node_filtered_networks$fsAndDiff_net,
  op_net,
  bouhaddou_randomRewiring_overlap$fsAndDiff_net,
  bouhaddou_randomRewiring_overlap$fsAndDiff_net_pval,
  title = 'Bouhaddou et al. (Random Rewiring)'
)

p <- ggarrange(wilkes_plot_randFiltering, wilkes_plot_randRewiring, bouhaddou_plot_randFiltering, bouhaddou_plot_randRewiring, 
          labels = c('A', '', 'B', ''),
          ncol = 1, nrow = 4, font.label=list(color="black",size=8))

ggsave('overlap_plot.svg', p, device = 'svg', units = 'mm', width = 55, height = 100)


## Plotting heatmap
#wilkes_randomFiltering_emp_p_vals <- c(
#  wilkes_randomFiltering_overlap[['diff_net_pval']],
#  wilkes_randomFiltering_overlap[['fs_net_pval']],
#  wilkes_randomFiltering_overlap[['fsAndDiff_net_pval']]
#) 
#wilkes_randomRewiring_emp_p_vals <- c(
#  wilkes_randomRewiring_overlap[['diff_net_pval']],
#  wilkes_randomRewiring_overlap[['fs_net_pval']],
#  wilkes_randomRewiring_overlap[['fsAndDiff_net_pval']]
#)

#bouhaddou_randomFiltering_emp_p_vals <- c(
#  bouhaddou_randomFiltering_overlap[['diff_net_pval']],
#  bouhaddou_randomFiltering_overlap[['fs_net_pval']],
#  bouhaddou_randomFiltering_overlap[['fsAndDiff_net_pval']]
#)
#bouhaddou_randomRewiring_emp_p_vals <- c(
#  bouhaddou_randomRewiring_overlap[['diff_net_pval']],
#  bouhaddou_randomRewiring_overlap[['fs_net_pval']],
#  bouhaddou_randomRewiring_overlap[['fsAndDiff_net_pval']]
#

#emp_p_vals <- data.frame(
#  'Wilkes et al. 2015' = wilkes_randomFiltering_emp_p_vals,
#  'Bouhaddou et al. 2023' = bouhaddou_randomFiltering_emp_p_vals,
#  'Wilkes et al. 2015' = wilkes_randomRewiring_emp_p_vals,
#  'Bouhaddou et al. 2023' = bouhaddou_randomRewiring_emp_p_vals,
#  check.names = F
#)
#rownames(emp_p_vals) <- c('DIFF', 'FS', 'DIFF + FS')

#cell_annotation <- data.frame(
#  rf1 = c('0.08', '< 0.001', '< 0.001'), 
#  rf2 = c('0.999', '< 0.001', '< 0.001'), 
#  rr1 = c('< 0.001', '<  0.001', '<  0.001'),
#  rr1 = c('< 0.001', '<  0.001', '<  0.001')
#)

#svg("NodeFilteredNetworks_pVals.svg")
#ComplexHeatmap::pheatmap(
#  emp_p_vals, 
#  color = colorRampPalette(c("navy", "white", "red"))(300),
#  breaks = seq(0, 0.1, length.out = 301),
#  cellwidth = 50,
#  cellheight = 50,
#  cluster_rows = F,
#  cluster_cols = F,
#  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#    grid.text(cell_annotation[i, j], x, y, gp = gpar(col = 'white'))
#  },
#  heatmap_legend_param = list(
#    title = 'P-value', at = c(0, 0.025, 0.05, 0.075, 0.1)
#  ),
#  column_split = c('Rand. Filtering', 'Rand. Filtering', 'Rand. Rewiring', 'Rand. Rewiring'),
#)
#dev.off()

## Wilkes PCST
wilkes_PCST_net <- read_tsv('./results/Wilkes2015/PCST_networks/MCF7-G2_vs_Mock_PCST_network.tsv', show_col_types = F)
wilkes_PCST_net$Interaction <- paste0(wilkes_PCST_net$Source, '_', wilkes_PCST_net$Target)

wilkes_PCST_FSfiltered_net <- read_tsv('./results/Wilkes2015/PCST_networks/MCF7-G2_vs_Mock_withFSfilter_PCST_network.tsv', show_col_types = F)
wilkes_PCST_FSfiltered_net$Interaction <- paste0(wilkes_PCST_FSfiltered_net$Source, '_', wilkes_PCST_FSfiltered_net$Target)

wilkes_PCST_Difffiltered_net <- read_tsv('./results/Wilkes2015/PCST_networks/MCF7-G2_vs_Mock_withDIFFfilter_PCST_network.tsv', show_col_types = F)
wilkes_PCST_Difffiltered_net$Interaction <- paste0(wilkes_PCST_Difffiltered_net$Source, '_', wilkes_PCST_Difffiltered_net$Target)

wilkes_PCST_FSandDiffFiltered_net <- read_tsv('./results/Wilkes2015/PCST_networks/MCF7-G2_vs_Mock_withFSandDIFFfilter_PCST_network.tsv', show_col_types = F)
wilkes_PCST_FSandDiffFiltered_net$Interaction <- paste0(wilkes_PCST_FSandDiffFiltered_net$Source, '_', wilkes_PCST_FSandDiffFiltered_net$Target)

wilkes_PCST_networks <- list(
  'PCST' = wilkes_PCST_net,
  'FS + PCST' = wilkes_PCST_FSfiltered_net,
  'DIFF + PCST' = wilkes_PCST_Difffiltered_net,
  'FS + DIFF + PCST' = wilkes_PCST_FSandDiffFiltered_net
)

wilkes_PCST_networks_randomFiltering_results <- list()
wilkes_PCST_networks_randomRewiring_results <- list()
for (filtered_net in names(wilkes_PCST_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    wilkes_sub_scores, 
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
bouhaddou_PCST_net <- read_tsv('./results/Bouhaddou2023/PCST_networks/VIC_10h_vs_Mock_10h_PCST_network.tsv', show_col_types = F)
bouhaddou_PCST_net$Interaction <- paste0(bouhaddou_PCST_net$Source, '_', bouhaddou_PCST_net$Target)

bouhaddou_PCST_FSfiltered_net <- read_tsv('./results/Bouhaddou2023/PCST_networks/VIC_10h_vs_Mock_10h_withFSfilter_PCST_network.tsv', show_col_types = F)
bouhaddou_PCST_FSfiltered_net$Interaction <- paste0(bouhaddou_PCST_FSfiltered_net$Source, '_', bouhaddou_PCST_FSfiltered_net$Target)

bouhaddou_PCST_Difffiltered_net <- read_tsv('./results/Bouhaddou2023/PCST_networks/VIC_10h_vs_Mock_10h_withDIFFfilter_PCST_network.tsv', show_col_types = F)
bouhaddou_PCST_Difffiltered_net$Interaction <- paste0(bouhaddou_PCST_Difffiltered_net$Source, '_', bouhaddou_PCST_Difffiltered_net$Target)

bouhaddou_PCST_FSandDiffFiltered_net <- read_tsv('./results/Bouhaddou2023/PCST_networks/VIC_10h_vs_Mock_10h_withFSandDIFFfilter_PCST_network.tsv', show_col_types = F)
bouhaddou_PCST_FSandDiffFiltered_net$Interaction <- paste0(bouhaddou_PCST_FSandDiffFiltered_net$Source, '_', bouhaddou_PCST_FSandDiffFiltered_net$Target)

bouhaddou_PCST_networks <- list(
  'PCST' = bouhaddou_PCST_net,
  'FS + PCST' = bouhaddou_PCST_FSfiltered_net,
  'DIFF + PCST' = bouhaddou_PCST_Difffiltered_net,
  'FS + DIFF + PCST' = bouhaddou_PCST_FSandDiffFiltered_net
)

bouhaddou_PCST_networks_randomFiltering_results <- list()
bouhaddou_PCST_networks_randomRewiring_results <- list()
for (filtered_net in names(bouhaddou_PCST_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    bouhaddou_sub_scores, 
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


bouhaddou_correlation_net <- read_tsv('./results/Bouhaddou2023/correlation_networks/VIC_10h.tsv', show_col_types = F)
# Only checking against E_data
bouhaddou_correlation_net$Source_LeadingProtein <- sapply(strsplit(bouhaddou_correlation_net$Source, '_'), function(x) { x[[1]] })
bouhaddou_correlation_net$Source <- bouhaddou_correlation_net$Source_LeadingProtein
bouhaddou_correlation_net$Interaction <- paste0(bouhaddou_correlation_net$Source, '_', bouhaddou_correlation_net$Target)

bouhaddou_node_filtered_networks$fs_net$Interaction <- paste0(bouhaddou_node_filtered_networks$fs_net$Source, '_', bouhaddou_node_filtered_networks$fs_net$Target)
bouhaddou_node_filtered_networks$diff_net$Interaction <- paste0(bouhaddou_node_filtered_networks$diff_net$Source, '_', bouhaddou_node_filtered_networks$diff_net$Target)
bouhaddou_node_filtered_networks$fsAndDiff_net$Interaction <- paste0(bouhaddou_node_filtered_networks$fsAndDiff_net$Source, '_', bouhaddou_node_filtered_networks$fsAndDiff_net$Target)

bouhaddou_corrFiltered_networks <- list(
  'corr_net' = unique(bouhaddou_correlation_net[, c('Source', 'Target')]),
  'corrAndDiffAndFs_net' = bouhaddou_node_filtered_networks$fsAndDiff_net[
    which(
      bouhaddou_node_filtered_networks$fsAndDiff_net$Interaction %in% bouhaddou_correlation_net$Interaction | 
      !(bouhaddou_node_filtered_networks$fsAndDiff_net$Source %in% bouhaddou_node_filtered_networks$fsAndDiff_net$Target_LeadingProtein)
    ), c('Source', 'Target')
  ]
)


## Bouhaddou CORR
bouhaddou_corrFilteredNets_randomFiltering_results <- list()
bouhaddou_corrFilteredNets_randomRewiring_results <- list()
for (filtered_net in names(bouhaddou_corrFiltered_networks)) {
  result_randomFiltering <- perform_randomFiltering_test(
    bouhaddou_sub_scores, 
    bouhaddou_corrFiltered_networks[[filtered_net]], 
    op_net
  )
  bouhaddou_corrFilteredNets_randomFiltering_results[[filtered_net]] <- result_randomFiltering[[1]]
  bouhaddou_corrFilteredNets_randomFiltering_results[[paste0(filtered_net, '_pval')]] <- result_randomFiltering[[2]]
  
  result_randomRewiring <- perform_randomRewiring_test(
    bouhaddou_corrFiltered_networks[[filtered_net]],
    op_net
  )
  bouhaddou_corrFilteredNets_randomRewiring_results[[filtered_net]] <- result_randomRewiring[[1]]
  bouhaddou_corrFilteredNets_randomRewiring_results[[paste0(filtered_net, '_pval')]] <- result_randomRewiring[[2]]
}


## plotting
wilkes_edgeFiltering_randomFiltering_emp_p_vals <- c(
  wilkes_PCST_networks_randomFiltering_results[['PCST_pval']],
  #wilkes_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  #wilkes_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  wilkes_PCST_networks_randomFiltering_results[['FS + DIFF + PCST_pval']],
  NA,
  #NA,
 # NA,
  NA
) 
wilkes_edgeFiltering_randomRewiring_emp_p_vals <- c(
  wilkes_PCST_networks_randomRewiring_results[['PCST_pval']],
  #wilkes_PCST_networks_randomRewiring_results[['FS + PCST_pval']],
  #wilkes_PCST_networks_randomRewiring_results[['DIFF + PCST_pval']],
  wilkes_PCST_networks_randomRewiring_results[['FS + DIFF + PCST_pval']],
  NA,
  #NA,
 # NA,
  NA
)

bouhaddou_edgeFiltering_randomFiltering_emp_p_vals <- c(
  bouhaddou_PCST_networks_randomFiltering_results[['PCST_pval']],
 # bouhaddou_PCST_networks_randomFiltering_results[['FS + PCST_pval']],
  #bouhaddou_PCST_networks_randomFiltering_results[['DIFF + PCST_pval']],
  bouhaddou_PCST_networks_randomFiltering_results[['FS + DIFF + PCST_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['corr_net_pval']],
#  bouhaddou_corrFilteredNets_randomFiltering_results[['corrAndFS_net_pval']],
  #bouhaddou_corrFilteredNets_randomFiltering_results[['corrAndDiff_net_pval']],
  bouhaddou_corrFilteredNets_randomFiltering_results[['corrAndDiffAndFs_net_pval']]
)

bouhaddou_edgeFiltering_randomRewiring_emp_p_vals <- c(
  bouhaddou_PCST_networks_randomRewiring_results[['PCST_pval']],
  #bouhaddou_PCST_networks_randomRewiring_results[['FS + PCST_pval']],
  #ouhaddou_PCST_networks_randomRewiring_results[['DIFF + PCST_pval']],
  bouhaddou_PCST_networks_randomRewiring_results[['FS + DIFF + PCST_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['corr_net_pval']],
#  bouhaddou_corrFilteredNets_randomRewiring_results[['corrAndFS_net_pval']],
  #bouhaddou_corrFilteredNets_randomRewiring_results[['corrAndDiff_net_pval']],
  bouhaddou_corrFilteredNets_randomRewiring_results[['corrAndDiffAndFs_net_pval']]
)

edgeFiltering_emp_p_vals <- data.frame(
  'Wilkes et al. 2015' = wilkes_edgeFiltering_randomFiltering_emp_p_vals,
  'Bouhaddou et al. 2023' = bouhaddou_edgeFiltering_randomFiltering_emp_p_vals,
  'Wilkes et al. 2015' = wilkes_edgeFiltering_randomRewiring_emp_p_vals,
  'Bouhaddou et al. 2023' = bouhaddou_edgeFiltering_randomRewiring_emp_p_vals,
  check.names = F
)

rownames(edgeFiltering_emp_p_vals) <- c('PCST', 'FS + DIFF + PCST', 'CORR', 'FS + DIFF + CORR')

cell_annotation <- data.frame(
  rf1 = c('0.07', '0.001', 'NA', 'NA'), 
  rf2 = c('0.005', '< 0.001', '< 0.001', '< 0.001'),
  rr1 = c('< 0.001', '< 0.001', 'NA', 'NA'),
  rr2 = c('< 0.001', '< 0.001', '< 0.001', '< 0.001')
)

svg("EdgeFilteredNetworks_pVals.svg")
ComplexHeatmap::pheatmap(
  edgeFiltering_emp_p_vals, 
  color = colorRampPalette(c("navy", "white", "red"))(500),
  breaks = seq(0, 0.1, length.out = 501),
  cellwidth = 60,
  cellheight = 60,
  cluster_rows = F,
  cluster_cols = F,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    text_color <- ifelse(is.na(edgeFiltering_emp_p_vals[i, j]), "black", "white")
    grid.text(cell_annotation[i, j], x, y, gp = gpar(col = text_color))
  },
  heatmap_legend_param = list(
    title = 'P-value', at = c(0, 0.025, 0.05, 0.075, 0.1)
  ),
  column_split = c('Random Filtering', 'Random Filtering', 'Random Rewiring', 'Random Rewiring'),
  row_split = c(rep('PCST', 2), rep('CORR', 2))
)
dev.off()


kinase_data <- load_kinase_data()



