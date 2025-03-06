library(tidyverse)
library(ComplexHeatmap)
library(igraph)
library(ggpubr)
library(gridExtra)
library(grid)
library(GGally)
library(patchwork)
# source('./src/UniprotIDMapping.R')

# Performs random filtering
perform_randomFiltering_test <- function(baseline_net, filtered_net, op_net, n_random = 1000) {
  baseline_net$Interaction <- paste0(baseline_net$Source, '_', baseline_net$Target)
  filtered_net$Interaction <- paste0(filtered_net$Source, '_', filtered_net$Target)
  filtered_net_interactions_in_op <- filtered_net[which(filtered_net$Interaction %in% op_net$Interaction),]
  no_randomInteractions <- c()
  for (i in 1:n_random) {
    net_random <- baseline_net[sample(nrow(baseline_net), nrow(filtered_net)), ]
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

# Performs random rewiring
perform_randomRewiring_test <- function(filtered_net, op_net, n_random = 1000) {
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

# Computes the randomization test results
get_filtered_KIN_results <- function(path, tested_alpha, tested_n, tested_beta, tested_gamma, tested_delta) {
  baseline_KINs <- list.files(paste0(path, 'baseline_KIN/'))
  baseline_KINs <- baseline_KINs[!grepl('baseline', baseline_KINs)]
  node_filtered_KINs <- list.files(paste0(path, 'node_filtered_KIN/'))
  edge_filtered_KINs <- list.files(paste0(path, 'edge_filtered_KIN/'))
  
  node_KINs <- list()
  pcst_KINs <- list()
  corr_x0_KINs <- list()
  corr_x1_KINs <- list()
  
  node_rf_results <- list()
  node_rr_results <- list()
  
  pcst_rr_results <- list()
  pcst_rf_results <- list()
  
  corr_x0_rr_results <- list()
  corr_x0_rf_results <- list()
  corr_x1_rr_results <- list()
  corr_x1_rf_results <- list()
  
  
  for (alpha in tested_alpha) {
    for (n in tested_n) {
      # get the baseline KIN
      baseline_KIN <- baseline_KINs[grepl(paste0('alpha', alpha, '_'), baseline_KINs)]
      baseline_KIN <- baseline_KIN[grepl(paste0('n', n, '_'), baseline_KIN)]
      baseline_KIN <- read_tsv(paste0(path, 'baseline_KIN/', baseline_KIN))
      baseline_KIN <- baseline_KIN[which(baseline_KIN$Edge_type == 'KS'),]
      
      for (beta in tested_beta) {
        for (gamma in tested_gamma) {
          # get the specific filtered KIN
          filtered_KIN <- node_filtered_KINs[grepl(paste0('alpha', alpha, '_'), node_filtered_KINs)] 
          filtered_KIN <- filtered_KIN[grepl(paste0('n', n, '_'), filtered_KIN)] 
          filtered_KIN <- filtered_KIN[grepl(paste0('gamma', gamma, '_'), filtered_KIN)]
          filtered_KIN <- filtered_KIN[grepl(paste0('beta', beta, '_'), filtered_KIN)] 
          filtered_KIN <- read_tsv(paste0(path, 'node_filtered_KIN/', filtered_KIN))
          
          # compute random filtering and random rewiring
          filtered_KIN <- filtered_KIN[which(filtered_KIN$Edge_type == 'KS'), ]
          node_KINs[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma)]] <- filtered_KIN
          #node_rf_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma)]] <- perform_randomFiltering_test(
          #  baseline_KIN, 
          #  filtered_KIN, 
          #  op_net
          #)
          #node_rr_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma)]] <- perform_randomRewiring_test(
          #  filtered_KIN, 
          #  op_net
          #)
          
          # get edge filtered KIN
          filtered_KIN <- edge_filtered_KINs[grepl(paste0('alpha', alpha, '_'), edge_filtered_KINs)]
          filtered_KIN <- filtered_KIN[grepl(paste0('n', n, '_'), filtered_KIN)] 
          filtered_KIN <- filtered_KIN[grepl(paste0('gamma', gamma, '_'), filtered_KIN)]
          filtered_KIN <- filtered_KIN[grepl(paste0('beta', beta, '_'), filtered_KIN)] 
        
          # get PCST filtered KIN
          filtered_PCST_KIN <- filtered_KIN[grepl('PCST_', filtered_KIN)]
          filtered_PCST_KIN <- read_tsv(paste0(path, 'edge_filtered_KIN/', filtered_PCST_KIN))
          filtered_PCST_KIN <- filtered_PCST_KIN[which(filtered_PCST_KIN$Edge_type == 'KS'), ]
          pcst_KINs[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma)]] <- filtered_PCST_KIN
          #pcst_rf_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma)]] <- perform_randomFiltering_test(
          #  baseline_KIN, 
          #  filtered_PCST_KIN, 
          #  op_net
          #)
          #pcst_rr_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma)]] <- perform_randomRewiring_test(
          #  filtered_PCST_KIN, 
          #  op_net
          #)
          
          # get CORR filtered KINs
          filtered_CORR_KIN <- filtered_KIN[grepl('CORR_', filtered_KIN)]
          if (length(filtered_CORR_KIN) > 0) {
            for (delta in tested_delta) {
              filtered_delta_CORR_KIN <- filtered_CORR_KIN[grepl(paste0('delta', delta, '_'), filtered_CORR_KIN)]
              filtered_CORR_x0_KIN <- filtered_delta_CORR_KIN[grepl('x0', filtered_delta_CORR_KIN)]
              filtered_CORR_x0_KIN <- read_tsv(paste0(path, 'edge_filtered_KIN/', filtered_CORR_x0_KIN))
              filtered_CORR_x0_KIN <- filtered_CORR_x0_KIN[which(filtered_CORR_x0_KIN$Edge_type == 'KS'),]
              corr_x0_KINs[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma, '_delta', delta)]] <- filtered_CORR_x0_KIN
              corr_x0_rr_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma, '_delta', delta)]] <- perform_randomFiltering_test(
                baseline_KIN, 
                filtered_CORR_x0_KIN, 
                op_net
              )
              corr_x0_rf_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma, '_delta', delta)]] <- perform_randomRewiring_test(
                filtered_CORR_x0_KIN, 
                op_net
              )
              
              filtered_CORR_x1_KIN <- filtered_delta_CORR_KIN[grepl('x1', filtered_delta_CORR_KIN)]
              filtered_CORR_x1_KIN <- read_tsv(paste0(path, 'edge_filtered_KIN/', filtered_CORR_x1_KIN))
              filtered_CORR_x1_KIN <- filtered_CORR_x1_KIN[which(filtered_CORR_x1_KIN$Edge_type == 'KS'),]
              corr_x1_KINs[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma, '_delta', delta)]] <- filtered_CORR_x1_KIN
              corr_x1_rr_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma, '_delta', delta)]] <- perform_randomFiltering_test(
                baseline_KIN, 
                filtered_CORR_x1_KIN, 
                op_net
              )
              corr_x1_rf_results[[paste0('alpha', alpha, '_n', n, '_beta', beta, '_gamma', gamma, '_delta', delta)]] <- perform_randomRewiring_test(
                filtered_CORR_x1_KIN, 
                op_net
              )
            }
          }
        }
      }
    }
  }
  results <- list(
    'node_KINs'=node_KINs,
    'pcst_KINs'=pcst_KINs,
    'corr_x0_KINs'=corr_x0_KINs,
    'corr_x1_KINs'=corr_x1_KINs,
    'node_rf'=node_rf_results,
    'node_rr'=node_rr_results,
    'pcst_rf'=pcst_rf_results,
    'pcst_rr'=pcst_rr_results,
    'corr_x0_rf'=corr_x0_rf_results,
    'corr_x0_rr'=corr_x0_rr_results,
    'corr_x1_rf'=corr_x1_rf_results,
    'corr_x1_rf'=corr_x1_rr_results
  )
  saveRDS(results, paste0(path, 'randomization_test_results.rds'))
  return(results)
}


get_plots <- function(results, tested_alpha, tested_n, tested_beta, tested_gamma, tested_delta, dataset='Wilkes2015') {

  plot_heatmap <- function(p_vals, KIN_rows, title) {
    p_vals_matrix <- matrix(unlist(p_vals), nrow = length(tested_alpha), byrow = TRUE)
    rownames(p_vals_matrix) <- paste0('alpha = ', tested_alpha)
    colnames(p_vals_matrix) <- sapply(strsplit(names(p_vals), '_'), function(x) paste0(x[[3]], ' ', x[[4]]))[1:(length(tested_n)*length(tested_beta)*length(tested_gamma))]
    
    p_vals_matrix[is.na(p_vals_matrix)] <- 1
    KIN_rows_matrix <- t(matrix(unlist(KIN_rows), nrow = length(tested_alpha), byrow = TRUE))
    
    heatmap <- Heatmap(t(p_vals_matrix), name = "p-values", 
              col = circlize::colorRamp2(c(0, 0.05, 0.1, 1), RColorBrewer::brewer.pal(4, 'YlOrRd')),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = rep(paste0('n = ', tested_n), each = ncol(p_vals_matrix)/length(tested_n)),
              cell_fun = function(j, i, x, y, width, height, fill) {
                value <- t(p_vals_matrix)[i, j]
                if (value == min(p_vals_matrix)) {
                  grid.text(paste0(sprintf("%.2e", value), ' (', KIN_rows_matrix[i, j], ')'), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
                } else {
                 grid.text(paste0(sprintf("%.2e", value), ' (', KIN_rows_matrix[i, j], ')'), x, y, gp = gpar(fontsize = 8))
                }

              },
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 45
    )
    
    heatmap <- draw(heatmap, 
         heatmap_legend_side = "right", 
         annotation_legend_side = "right",
         column_title = title)
    
    return(heatmap)  
  }
  
  heatmaps <- list()
  
  # node rr filtering 
  rf_p_vals <- lapply(results$node_rf, function(x) x$p_val)
  rr_p_vals <- lapply(results$node_rf, function(x) x$p_val)
  KIN_rows <- lapply(results$node_KINs, function(x) nrow(x))
  heatmaps[['rf_node_filtering']] <- plot_heatmap(rf_p_vals, KIN_rows, title = ifelse(dataset=='Wilkes2015', 'Random filtering: PI3K inhibition in resistant breast cancer cell line (DIFF + FS)', 'Random filtering: SARS-CoV-2 (DIFF + FS)'))
  heatmaps[['rr_node_filtering']] <- plot_heatmap(rr_p_vals, KIN_rows, title = ifelse(dataset=='Wilkes2015', 'Random rewiring: I3K inhibition in resistant breast cancer cell line (DIFF + FS)', 'Random rewiring: SARS-CoV-2 (DIFF + FS)'))
  # pcst filtering
  rf_p_vals <- lapply(results$pcst_rf, function(x) x$p_val)
  rr_p_vals <- lapply(results$pcst_rf, function(x) x$p_val)
  KIN_rows <- lapply(results$pcst_KINs, function(x) nrow(x))
  heatmaps[['rf_pcst_filtering']] <- plot_heatmap(rf_p_vals, KIN_rows, title = ifelse(dataset=='Wilkes2015', 'Random filtering: PI3K inhibition in resistant breast cancer cell line (DIFF + FS + PCST)', 'Random filtering: SARS-CoV-2 (DIFF + FS + PCST)'))
  heatmaps[['rr_pcst_filtering']] <- plot_heatmap(rr_p_vals, KIN_rows, title = ifelse(dataset=='Wilkes2015', 'Random rewiring: PI3K inhibition in resistant breast cancer cell line (DIFF + FS + PCST)', 'Random rewiring: SARS-CoV-2 (DIFF + FS + PCST)'))
  
  return(heatmaps)
}



# Loading Omnipath network
op_net <- read_csv('../data/OmniPath/2023_10_11_KinaseDataOmniPath.csv', show_col_types = F)
op_net$Interaction <- paste0(op_net$enzyme, '_', op_net$TARGET_UP_ID, '_', op_net$TARGET_RES, op_net$TARGET_POS)

# Tested hyperparameters
tested_alpha <- c(0.85, 0.9, 0.95)
tested_n <- c(10, 15, 20)
tested_beta <- c(0.0, 0.2, 0.4, 0.6)
tested_gamma <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5)
tested_delta <- c(0.7, 0.8, 0.9)

# Analysis for Wilkes et al. 2015
wilkes_path <- '../results/Hyperparameters/Wilkes2015/'
wilkes_filtered_KIN_results <- get_filtered_KIN_results(
  path = wilkes_path, 
  tested_alpha,
  tested_n,
  tested_beta,
  tested_gamma,
  tested_delta
)

wilkes_filtered_KIN_results <- readRDS(paste0(wilkes_path, 'randomization_test_results.rds'))
wilkes_plots <- get_plots(wilkes_filtered_KIN_results, tested_alpha, tested_n, tested_beta, tested_gamma, tested_delta, dataset='Wilkes2015') 
# Analysis for Bouhaddou et al. 2023
bouhaddou_path <- '../results/Hyperparameters/Bouhaddou2023/'
bouhaddou_filtered_KIN_results <- get_filtered_KIN_results(
   path = bouhaddou_path,
   tested_alpha,
   tested_n,
   tested_beta,
   tested_gamma,
   tested_delta
)
# bouhaddou_filtered_KIN_results <- readRDS(paste0(bouhaddou_path, 'randomization_test_results.rds'))
# bouhaddou_plots <- get_plots(bouhaddou_filtered_KIN_results, tested_alpha, tested_n, tested_beta, tested_gamma, tested_delta, dataset='Bouhaddou2023')



plot_hyperparameter_effect_line <- function(wilkes_p_vals, bouhaddou_p_vals, hyperparameter, hyperparameter_values, y_lab="Mean P-value", legend="") {
  names(wilkes_p_vals) <- paste0(names(wilkes_p_vals), '_')
  names(bouhaddou_p_vals) <- paste0(names(bouhaddou_p_vals), '_')
  print(names(bouhaddou_p_vals))
  data <- data.frame()
  for (param in hyperparameter_values) {
    print(paste0(hyperparameter,param, '_'))
    p_vals <- unlist(wilkes_p_vals[names(wilkes_p_vals)[grepl(paste0(hyperparameter,param, '_'), names(wilkes_p_vals))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- paste('PI3K inhibition in resistant breast cancer cell line', legend) 
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition))
  }  
  for (param in hyperparameter_values) {
    p_vals <- unlist(bouhaddou_p_vals[names(bouhaddou_p_vals)[grepl(paste0(hyperparameter,param, '_'), names(bouhaddou_p_vals))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- paste('SARS-CoV-2', legend) 
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition))
  }  
  print(data)
  data$hyperparam <- as.factor(data$hyperparam)
  if (hyperparameter == 'alpha') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      bquote(alpha == .(x))  # Uses bquote for dynamic labels
    })
  } else if (hyperparameter == 'n') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      paste0('n = ', x)  # Uses bquote for dynamic labels
    })
  } else if (hyperparameter == 'beta') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      bquote(beta == .(x))  # Uses bquote for dynamic labels
    })
  } else if (hyperparameter == 'gamma') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      bquote(gamma == .(x))  # Uses bquote for dynamic labels
    }) 
  }
  if (hyperparameter == 'n') {
    plot <- ggplot(data, aes(x = hyperparam, y = mean_p, group = condition, color = condition)) +
      geom_line(aes(group=condition)) +  # Line plot for mean p-value
      geom_point(aes(group=condition)) +  # Points for means
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = condition), alpha = 0.2, color = NA) +  # Confidence interval
      scale_x_discrete(labels = hyperparam_labels) +
      # scale_y_continuous(limits = c(0, round(max(data$upper) + 0.05, str_count(gsub(".*[.]","",as.character(max(data$upper))), "0")+1)), breaks = waiver(), n.breaks = 5) +
      xlab("Hyperparameter values") +
      ylab(y_lab) + 
      theme_minimal()
  } else {
    plot <- ggplot(data, aes(x = hyperparam, y = mean_p, group = condition, color = condition)) +
      geom_line(aes(group=condition)) +  # Line plot for mean p-value
      geom_point(aes(group=condition)) +  # Points for means
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = condition), alpha = 0.2, color = NA) +  # Confidence interval
      scale_x_discrete(labels = do.call(expression, hyperparam_labels)) +
      xlab("Hyperparameter values") +
      ylab(y_lab) + 
      theme_minimal()
  }
  return(plot)
}

bouhaddou_filtered_KIN_results <- readRDS(paste0(bouhaddou_path, 'randomization_test_results.rds'))
wilkes_filtered_KIN_results <- readRDS(paste0(wilkes_path, 'randomization_test_results.rds'))

# NODE RF plots
wilkes_node_rf_p_vals <- lapply(wilkes_filtered_KIN_results$node_rf, function(x) x$p_val)
bouhaddou_node_rf_p_vals <- lapply(bouhaddou_filtered_KIN_results$node_rf, function(x) x$p_val)
alpha_node_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, "alpha", tested_alpha)
n_node_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, 'n', tested_n)
gamma_node_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, 'gamma', tested_gamma)
beta_node_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, 'beta', tested_beta)

node_rf_plot <- alpha_node_rf_effect_plot + 
  n_node_rf_effect_plot + 
  gamma_node_rf_effect_plot + 
  beta_node_rf_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Node RF (DIFF, FS)") & 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

# NODE RR plots
wilkes_node_rr_p_vals <- lapply(wilkes_filtered_KIN_results$node_rr, function(x) x$p_val) 
bouhaddou_node_rr_p_vals <- lapply(bouhaddou_filtered_KIN_results$node_rr, function(x) x$p_val) 
alpha_node_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, "alpha", tested_alpha)
n_node_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, 'n', tested_n, y_lab="")
gamma_node_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, 'gamma', tested_gamma)
beta_node_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, 'beta', tested_beta, y_lab="")

node_rr_plot <- alpha_node_rr_effect_plot + 
  n_node_rr_effect_plot + 
  gamma_node_rr_effect_plot + 
  beta_node_rr_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Node RR (DIFF, FS)") & 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

# PCST RF plots
wilkes_pcst_rf_p_vals <- lapply(wilkes_filtered_KIN_results$pcst_rf, function(x) x$p_val)
bouhaddou_pcst_rf_p_vals <- lapply(bouhaddou_filtered_KIN_results$pcst_rf, function(x) x$p_val)
alpha_pcst_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, "alpha", tested_alpha)
n_pcst_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, 'n', tested_n, y_lab="")
gamma_pcst_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, 'gamma', tested_gamma)
beta_pcst_rf_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, 'beta', tested_beta, y_lab="")

pcst_rf_plot <- alpha_pcst_rf_effect_plot + 
  n_pcst_rf_effect_plot + 
  gamma_pcst_rf_effect_plot + 
  beta_pcst_rf_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Edge RF (PCST, DIFF, FS)") & 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

# PCST RR plots
wilkes_pcst_rr_p_vals <- lapply(wilkes_filtered_KIN_results$pcst_rr, function(x) x$p_val)
bouhaddou_pcst_rr_p_vals <- lapply(bouhaddou_filtered_KIN_results$pcst_rr, function(x) x$p_val)
alpha_pcst_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, "alpha", tested_alpha)
n_pcst_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, 'n', tested_n)
gamma_pcst_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, 'gamma', tested_gamma)
beta_pcst_rr_effect_plot <- plot_hyperparameter_effect_line(wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, 'beta', tested_beta)

pcst_rr_plot <- alpha_pcst_rr_effect_plot + 
  n_pcst_rr_effect_plot + 
  gamma_pcst_rr_effect_plot + 
  beta_pcst_rr_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Edge RR (PCST, DIFF, FS)") & 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))


# Effect size plots

# Node filters
wilkes_KIN_node_rows <- lapply(wilkes_filtered_KIN_results$node_KINs, function(x) nrow(x))
bouhaddou_KIN_node_rows <- lapply(bouhaddou_filtered_KIN_results$node_KINs, function(x) nrow(x))
alpha_node_edge_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, "alpha", tested_alpha, y_lab = "Number of edges")
n_node_edge_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, 'n', tested_n, y_lab = "Number of edges")
gamma_node_edge_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, 'gamma', tested_gamma, y_lab = "Number of edges")
beta_node_edge_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, 'beta', tested_beta, y_lab = "Number of edges")

node_edge_effect_plot <- alpha_node_edge_effect_plot + 
  n_node_edge_effect_plot + 
  gamma_node_edge_effect_plot + 
  beta_node_edge_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Number of edges  (DIFF, FS)") & 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))


wilkes_KIN_edge_rows <- lapply(wilkes_filtered_KIN_results$pcst_KINs, function(x) nrow(x))
bouhaddou_KIN_edge_rows <- lapply(bouhaddou_filtered_KIN_results$pcst_KINs, function(x) nrow(x))
alpha_pcst_edge_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, "alpha", tested_alpha, y_lab = "Number of edges", legend = "(with PCST)")
n_pcst_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, 'n', tested_n, y_lab = "Number of edges", legend = "(with PCST)")
gamma_pcst_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, 'gamma', tested_gamma, y_lab = "Number of edges", legend = "(with PCST)")
beta_pcst_effect_plot <- plot_hyperparameter_effect_line(wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, 'beta', tested_beta, y_lab = "Number of edges", legend = "(with PCST)")
pcst_edge_effect_plot <- alpha_pcst_edge_effect_plot + 
  n_pcst_effect_plot + 
  gamma_pcst_effect_plot + 
  beta_pcst_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Number of edges  (PCST, DIFF, FS)") & 
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))



plot_total_effect_size <- function(wilkes_p_vals_noderf, bouhaddou_p_vals_noderf, wilkes_p_vals_pcstrf, bouhaddou_p_vals_pcstrf, hyperparameter, hyperparameter_values, y_lab="Mean P-value") {
  names(wilkes_p_vals_noderf) <- paste0(names(wilkes_p_vals_noderf), '_')
  names(bouhaddou_p_vals_noderf) <- paste0(names(bouhaddou_p_vals_noderf), '_')
  names(wilkes_p_vals_pcstrf) <- paste0(names(wilkes_p_vals_pcstrf), '_')
  names(bouhaddou_p_vals_pcstrf) <- paste0(names(bouhaddou_p_vals_pcstrf), '_')
  data <- data.frame()
  for (param in hyperparameter_values) {
    p_vals <- unlist(wilkes_p_vals_noderf[names(wilkes_p_vals_noderf)[grepl(paste0(hyperparameter,param, '_'), names(wilkes_p_vals_noderf))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- 'PI3K inhibition in resistant breast cancer cell line (without PCST)' 
    col <- 'PI3K inhibition in resistant breast cancer cell line'
    linetype <- 'without PCST'
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition, col = col, linetype = linetype))
  }  
  for (param in hyperparameter_values) {
    p_vals <- unlist(bouhaddou_p_vals_noderf[names(bouhaddou_p_vals_noderf)[grepl(paste0(hyperparameter,param, '_'), names(bouhaddou_p_vals_noderf))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- 'SARS-CoV-2 (without PCST)'
    col <- 'SARS-CoV-2'
    linetype <- 'without PCST'
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition, col = col, linetype = linetype))
  }  
  for (param in hyperparameter_values) {
    p_vals <- unlist(wilkes_p_vals_pcstrf[names(wilkes_p_vals_pcstrf)[grepl(paste0(hyperparameter,param, '_'), names(wilkes_p_vals_pcstrf))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- 'PI3K inhibition in resistant breast cancer cell line (with PCST)' 
    col <- 'PI3K inhibition in resistant breast cancer cell line'
    linetype <- 'with PCST'
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition, col = col, linetype = linetype))
  }  
  for (param in hyperparameter_values) {
    p_vals <- unlist(bouhaddou_p_vals_pcstrf[names(bouhaddou_p_vals_pcstrf)[grepl(paste0(hyperparameter,param, '_'), names(bouhaddou_p_vals_pcstrf))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- 'SARS-CoV-2 (with PCST)'
    col <- 'SARS-CoV-2'
    linetype <- 'with PCST'
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition, col = col, linetype = linetype))
  } 
  print(data)
  data$hyperparam <- as.factor(data$hyperparam)
  if (hyperparameter == 'alpha') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      bquote(alpha == .(x))  # Uses bquote for dynamic labels
    })
  } else if (hyperparameter == 'n') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      paste0('n = ', x)  # Uses bquote for dynamic labels
    })
  } else if (hyperparameter == 'beta') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      bquote(beta == .(x))  # Uses bquote for dynamic labels
    })
  } else if (hyperparameter == 'gamma') {
    hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
      bquote(gamma == .(x))  # Uses bquote for dynamic labels
    }) 
  }
  if (hyperparameter == 'n') {
    plot <- ggplot(data, aes(x = hyperparam, y = mean_p, group = condition, color = col, linetype = linetype)) +
      geom_line(aes(group=condition)) +  # Line plot for mean p-value
      geom_point(aes(group=condition)) +  # Points for means
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = col), alpha = 0.2, color = NA) +  # Confidence interval
      scale_x_discrete(labels = hyperparam_labels) +
      # scale_y_continuous(limits = c(0, round(max(data$upper) + 0.05, str_count(gsub(".*[.]","",as.character(max(data$upper))), "0")+1)), breaks = waiver(), n.breaks = 5) +
      xlab("Hyperparameter values") +
      ylab(y_lab) + 
      theme_minimal()
  } else {
    plot <- ggplot(data, aes(x = hyperparam, y = mean_p, group = condition, color = col, linetype = linetype)) +
      geom_line(aes(group=condition)) +  # Line plot for mean p-value
      geom_point(aes(group=condition)) +  # Points for means
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = col), alpha = 0.2, color = NA) +  # Confidence interval
      scale_x_discrete(labels = do.call(expression, hyperparam_labels)) +
      xlab("Hyperparameter values") +
      ylab(y_lab) + 
      theme_minimal()
  }
  return(plot)
}

# EFFECT Size plot
wilkes_KIN_node_rows <- lapply(wilkes_filtered_KIN_results$node_KINs, function(x) nrow(x))
bouhaddou_KIN_node_rows <- lapply(bouhaddou_filtered_KIN_results$node_KINs, function(x) nrow(x))
wilkes_KIN_edge_rows <- lapply(wilkes_filtered_KIN_results$pcst_KINs, function(x) nrow(x))
bouhaddou_KIN_edge_rows <- lapply(bouhaddou_filtered_KIN_results$pcst_KINs, function(x) nrow(x))

alpha_node_edge_effect_plot <- plot_total_effect_size(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, "alpha", tested_alpha, y_lab = "Number of edges")
n_node_edge_effect_plot <- plot_total_effect_size(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, 'n', tested_n, y_lab = "Number of edges")
gamma_node_edge_effect_plot <- plot_total_effect_size(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, 'gamma', tested_gamma, y_lab = "Number of edges")
beta_node_edge_effect_plot <- plot_total_effect_size(wilkes_KIN_node_rows, bouhaddou_KIN_node_rows, wilkes_KIN_edge_rows, bouhaddou_KIN_edge_rows, 'beta', tested_beta, y_lab = "Number of edges")

node_edge_effect_plot <- alpha_node_edge_effect_plot + 
  n_node_edge_effect_plot + 
  gamma_node_edge_effect_plot + 
  beta_node_edge_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Number of edges") &
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm')) &
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))



## RF Plot
wilkes_node_rf_p_vals <- lapply(wilkes_filtered_KIN_results$node_rf, function(x) x$p_val)
bouhaddou_node_rf_p_vals <- lapply(bouhaddou_filtered_KIN_results$node_rf, function(x) x$p_val)
wilkes_pcst_rf_p_vals <- lapply(wilkes_filtered_KIN_results$pcst_rf, function(x) x$p_val)
bouhaddou_pcst_rf_p_vals <- lapply(bouhaddou_filtered_KIN_results$pcst_rf, function(x) x$p_val)

alpha_rf_effect_plot <- plot_total_effect_size(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, "alpha", tested_alpha)
n_rf_effect_plot <- plot_total_effect_size(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, 'n', tested_n, y_lab="")
gamma_rf_effect_plot <- plot_total_effect_size(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, 'gamma', tested_gamma)
beta_rf_effect_plot <- plot_total_effect_size(wilkes_node_rf_p_vals, bouhaddou_node_rf_p_vals, wilkes_pcst_rf_p_vals, bouhaddou_pcst_rf_p_vals, 'beta', tested_beta, y_lab="")

rf_test_plot <- alpha_rf_effect_plot + 
  n_rf_effect_plot + 
  gamma_rf_effect_plot + 
  beta_rf_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Random filtering test") &
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm')) &
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

manuscript_hyperparameter_plot <- grid.arrange(rf_test_plot, node_edge_effect_plot, nrow = 2)

## RR Plot
wilkes_node_rr_p_vals <- lapply(wilkes_filtered_KIN_results$node_rr, function(x) x$p_val)
bouhaddou_node_rr_p_vals <- lapply(bouhaddou_filtered_KIN_results$node_rr, function(x) x$p_val)
wilkes_pcst_rr_p_vals <- lapply(wilkes_filtered_KIN_results$pcst_rr, function(x) x$p_val)
bouhaddou_pcst_rr_p_vals <- lapply(bouhaddou_filtered_KIN_results$pcst_rr, function(x) x$p_val)

alpha_rr_effect_plot <- plot_total_effect_size(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, "alpha", tested_alpha)
n_rr_effect_plot <- plot_total_effect_size(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, 'n', tested_n, y_lab="")
gamma_rr_effect_plot <- plot_total_effect_size(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, 'gamma', tested_gamma)
beta_rr_effect_plot <- plot_total_effect_size(wilkes_node_rr_p_vals, bouhaddou_node_rr_p_vals, wilkes_pcst_rr_p_vals, bouhaddou_pcst_rr_p_vals, 'beta', tested_beta, y_lab="")

rr_test_plot <- alpha_rr_effect_plot + 
  n_rr_effect_plot + 
  gamma_rr_effect_plot + 
  beta_rr_effect_plot + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = "Hyperparameter evaluation: Random rewiring test") &
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm')) &
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))


## CORR plot
bouhaddou_corr_KIN_files <- list.files('./results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/')  
bouhaddou_corr_KIN_files <- bouhaddou_corr_KIN_files[grepl('x0', bouhaddou_corr_KIN_files)]
bouhaddou_corr_KIN_sizes <- list()
for (f in bouhaddou_corr_KIN_files) {
  d <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/', f))
  d <- d[d$Edge_type == 'KS',]
  bouhaddou_corr_KIN_sizes[[f]] <- nrow(d)
}

bouhaddou_node_KIN_files <- list.files('./results/Hyperparameters/Bouhaddou2023/node_filtered_KIN/')
bouhaddou_node_KIN_files <- bouhaddou_node_KIN_files[grepl('_KIN.tsv', bouhaddou_node_KIN_files)]
bouhaddou_node_KIN_sizes <- list()
for (f in bouhaddou_node_KIN_files) {
  d <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/node_filtered_KIN/', f))
  d <- d[d$Edge_type == 'KS',]
  bouhaddou_node_KIN_sizes[[f]] <- nrow(d)
}

corr_plot <- function(corr_KIN_sizes, node_KIN_sizes, hyperparameter='delta', hyperparameter_values=c(0.7,0.8,0.9)) {
  data <- data.frame()
  for (param in hyperparameter_values) {
    p_vals <- unlist(corr_KIN_sizes[names(corr_KIN_sizes)[grepl(paste0(hyperparameter,param, '_'), names(corr_KIN_sizes))]])
    p_vals <- p_vals[!is.na(p_vals)]
    mean_p <- mean(p_vals)
    se <- sd(p_vals) / sqrt(length(p_vals))
    lower <- max(c(mean_p - 1.96 * se, 0))
    upper <- mean_p + 1.96 * se
    condition <- 'SARS-CoV-2'
    col <- 'SARS-CoV-2'
    linetype <- ''
    data <- rbind(data, list(hyperparam = param, mean_p = mean_p, lower = lower, upper = upper, condition = condition, col = col, linetype = linetype))
  }
  node_KIN_sizes <- unlist(node_KIN_sizes)
  mean_without_delta <- mean(node_KIN_sizes)
  se_without_delta <- sd(node_KIN_sizes) / sqrt(length(node_KIN_sizes))
  lower <- max(c(mean_without_delta - 1.96 * se_without_delta, 0))
  upper <- mean_without_delta + 1.96 * se_without_delta
  
  data <- rbind(data, list(hyperparam = '0', mean_p = mean_without_delta, lower = lower, upper = upper, condition = 'SARS-CoV-2', col = 'SARS-CoV-2', linetype = ''))
  data$hyperparam <- as.factor(data$hyperparam)
  hyperparam_labels <- sapply(levels(data$hyperparam), function(x) {
    bquote(delta == .(x))  # Uses bquote for dynamic labels
  })
  plot <- ggplot(data, aes(x = hyperparam, y = mean_p, group = condition, color = col)) +
    geom_line(aes(group=condition)) +  # Line plot for mean p-value
    geom_point(aes(group=condition)) +  # Points for means
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = col), alpha = 0.2, color = NA) +  # Confidence interval
    scale_x_discrete(labels = do.call(expression, hyperparam_labels)) +
    scale_color_discrete(name = NULL) +  # Remove legend title for color
    scale_fill_discrete(name = NULL) +
    xlab("Hyperparameter values") +
    ylab("Number of edges") + 
    theme_minimal()
}
  
  



