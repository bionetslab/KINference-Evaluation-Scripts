library(tidyverse)
library(ggplot2)
library(readxl)

biomart <- read_tsv('../data/biomart/mart_export_human_uniprot.txt')

gene2pubmed <- read_tsv('../data/gene2pubmed/gene2pubmed.tsv')
gene2ensembl <- read_tsv('../data/gene2pubmed/gene2ensembl.tsv')
gene2pubmed <- gene2pubmed[gene2pubmed[['#tax_id']] == 9606, ]
gene2ensembl <- gene2ensembl[gene2ensembl[['#tax_id']] == 9606, ]

gene2pubmed_distribution <- table(gene2pubmed$GeneID)

df <- tibble(
  "Gene" = names(gene2pubmed_distribution),
  "Occurences" = as.integer(unname(gene2pubmed_distribution))
)

get_pubmed_distribution <- function(data) {
  targets <- unique(data$Target_Uniprot)
  targets <- targets[!is.na(targets)]
  ensg_IDs <- unique(biomart[which(biomart$`UniProtKB Gene Name ID` %in% targets), ][['Gene stable ID']])
  geneIDs <- unique(gene2ensembl[which(gene2ensembl$Ensembl_gene_identifier %in% ensg_IDs),]$'GeneID')
  distribution <- table(gene2pubmed[which(gene2pubmed$GeneID %in% geneIDs),]$GeneID)
  df <- tibble(
    "Gene" = names(distribution),
    "Occurences" = as.integer(unname(distribution))
  )
  return(df)
}



bouhaddou_baselineKIN <- read_tsv('../results/Hyperparameters/Bouhaddou2023/baseline_KIN/Bouhaddou2023_n15_alpha0.9_KIN.tsv')
bouhaddou_baselineKIN_df <- get_pubmed_distribution(bouhaddou_baselineKIN)

wilkes_baseline_KIN <- read_tsv('../results/Hyperparameters/Wilkes2015/baseline_KIN/Wilkes2015_n15_alpha0.9_KIN.tsv')
wilkes_baselineKIN_df <- get_pubmed_distribution(wilkes_baseline_KIN)

tested_gamma <- c(0.5, 1, 1.5, 2, 2.5)
tested_beta <- c(0.2, 0.4, 0.6)
tested_delta <- c(0.7, 0.8, 0.9)


## BOUHADDOU et al.
bouhaddou_node_literatureBias_gamma <- list()
bouhaddou_pcst_literatureBias_gamma <- list()
for (g in tested_gamma) {
  node_DKIN <- read_tsv(paste0('../results/Hyperparameters/Bouhaddou2023/node_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma', g, '_beta0_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('../results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma', g, '_beta0_PCST_filtered_KIN.tsv'))
  
  bouhaddou_node_literatureBias_gamma[[as.character(g)]] <- get_pubmed_distribution(node_DKIN)$Occurences
  bouhaddou_pcst_literatureBias_gamma[[as.character(g)]] <- get_pubmed_distribution(pcst_DKIN)$Occurences
}

bouhaddou_node_literatureBias_gamma_plot <- tibble(
  Values = c(bouhaddou_baselineKIN_df$Occurences, unname(unlist(bouhaddou_node_literatureBias_gamma))),
  Group = c(rep('Baseline KIN', nrow(bouhaddou_baselineKIN_df)), rep(names(bouhaddou_node_literatureBias_gamma), lengths(bouhaddou_node_literatureBias_gamma)))
)
bouhaddou_pcst_literatureBias_gamma_plot <- tibble(
  Values = c(bouhaddou_baselineKIN_df$Occurences, unname(unlist(bouhaddou_pcst_literatureBias_gamma))),
  Group = c(rep('Baseline KIN', nrow(bouhaddou_baselineKIN_df)), rep(names(bouhaddou_pcst_literatureBias_gamma), lengths(bouhaddou_pcst_literatureBias_gamma)))
)

# Convert Group labels to expressions for proper formatting
bouhaddou_node_literatureBias_gamma_plot$Group <- factor(
  bouhaddou_node_literatureBias_gamma_plot$Group,
  levels = c("Baseline KIN", names(bouhaddou_node_literatureBias_gamma)),
  labels = c("Baseline~KIN", paste0("gamma == ", names(bouhaddou_node_literatureBias_gamma)))
)
bouhaddou_pcst_literatureBias_gamma_plot$Group <- factor(
  bouhaddou_pcst_literatureBias_gamma_plot$Group,
  levels = c("Baseline KIN", names(bouhaddou_pcst_literatureBias_gamma)),
  labels = c("Baseline~KIN", paste0("gamma == ", names(bouhaddou_pcst_literatureBias_gamma)))
)

p1_gamma <- ggplot(bouhaddou_node_literatureBias_gamma_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Without PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  scale_fill_manual(
    values = c("steelblue", "tomato", "darkgreen", "purple", "orange")  # Adjust colors if needed
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

p2_gamma <- ggplot(bouhaddou_pcst_literatureBias_gamma_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "With PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  scale_fill_manual(
    values = c("steelblue", "tomato", "darkgreen", "purple", "orange")  # Adjust colors if needed
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

bouhaddou_node_literatureBias_beta <- list()
bouhaddou_pcst_literatureBias_beta <- list()
for (b in tested_beta) {
  node_DKIN <- read_tsv(paste0('../results/Hyperparameters/Bouhaddou2023/node_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma0_beta', b, '_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('../results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma0_beta', b, '_PCST_filtered_KIN.tsv'))
  
  bouhaddou_node_literatureBias_beta[[as.character(b)]] <- get_pubmed_distribution(node_DKIN)$Occurences
  bouhaddou_pcst_literatureBias_beta[[as.character(b)]] <- get_pubmed_distribution(pcst_DKIN)$Occurences
}

bouhaddou_node_literatureBias_beta_plot <- tibble(
  Values = c(bouhaddou_baselineKIN_df$Occurences, unname(unlist(bouhaddou_node_literatureBias_beta))),
  Group = c(rep('Baseline KIN', nrow(bouhaddou_baselineKIN_df)), rep(names(bouhaddou_node_literatureBias_beta), lengths(bouhaddou_node_literatureBias_beta)))
)
bouhaddou_pcst_literatureBias_beta_plot <- tibble(
  Values = c(bouhaddou_baselineKIN_df$Occurences, unname(unlist(bouhaddou_pcst_literatureBias_beta))),
  Group = c(rep('Baseline KIN', nrow(bouhaddou_baselineKIN_df)), rep(names(bouhaddou_pcst_literatureBias_beta), lengths(bouhaddou_pcst_literatureBias_beta)))
)

# Convert Group labels to expressions for proper formatting
bouhaddou_node_literatureBias_beta_plot$Group <- factor(
  bouhaddou_node_literatureBias_beta_plot$Group,
  levels = c("Baseline KIN", names(bouhaddou_node_literatureBias_beta)),
  labels = c("Baseline~KIN", paste0("beta == ", names(bouhaddou_node_literatureBias_beta)))
)
bouhaddou_pcst_literatureBias_beta_plot$Group <- factor(
  bouhaddou_pcst_literatureBias_beta_plot$Group,
  levels = c("Baseline KIN", names(bouhaddou_pcst_literatureBias_beta)),
  labels = c("Baseline~KIN", paste0("beta == ", names(bouhaddou_pcst_literatureBias_beta)))
)

p1_beta <- ggplot(bouhaddou_node_literatureBias_beta_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Without PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  #scale_fill_manual(
  #  values = c("steelblue", "tomato", "darkgreen", "purple")  # Adjust colors if needed
  #) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

p2_beta <- ggplot(bouhaddou_pcst_literatureBias_beta_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "With PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  #scale_fill_manual(
  #  values = c("steelblue", "tomato", "darkgreen", "purple")  # Adjust colors if needed
  #) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

bouhaddou_literatureBias_delta <- list()
for (d in tested_delta) {
  DKIN <- read_tsv(paste0('../results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n10_alpha0.85_gamma0_beta0_delta', d, '_CORR_filtered_x0_KIN.tsv'))
  bouhaddou_literatureBias_delta[[as.character(d)]] <- get_pubmed_distribution(DKIN)$Occurences
}

bouhaddou_literatureBias_delta_plot <- tibble(
  Values = c(bouhaddou_baselineKIN_df$Occurences, unname(unlist(bouhaddou_literatureBias_delta))),
  Group = c(rep('Baseline KIN', nrow(bouhaddou_baselineKIN_df)), rep(names(bouhaddou_literatureBias_delta), lengths(bouhaddou_literatureBias_delta)))
)

# Convert Group labels to expressions for proper formatting
bouhaddou_literatureBias_delta_plot$Group <- factor(
  bouhaddou_literatureBias_delta_plot$Group,
  levels = c("Baseline KIN", names(bouhaddou_literatureBias_delta)),
  labels = c("Baseline~KIN", paste0("delta == ", names(bouhaddou_literatureBias_delta)))
)

p1_delta <- ggplot(bouhaddou_literatureBias_delta_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = expression(paste("Hyperparameter evaluation: Literature bias of ", delta)), x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  scale_fill_manual(
    values = c("steelblue", "tomato", "darkgreen", "purple", "orange")  # Adjust colors if needed
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(labels = scales::parse_format())


gamma_plot <- p1_gamma + 
  p2_gamma + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = expression(paste("Hyperparameter evaluation: Literature bias of ", gamma))) &
  theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm'))

beta_plot <- p1_beta + 
  p2_beta + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = expression(paste("Hyperparameter evaluation: Literature bias of ", beta))) &
  theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm'))


## WILKES et al.
wilkes_node_literatureBias_gamma <- list()
wilkes_pcst_literatureBias_gamma <- list()
for (g in tested_gamma) {
  node_DKIN <- read_tsv(paste0('../results/Hyperparameters/Wilkes2015/node_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma', g, '_beta0_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('../results/Hyperparameters/Wilkes2015/edge_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma', g, '_beta0_PCST_filtered_KIN.tsv'))
  
  wilkes_node_literatureBias_gamma[[as.character(g)]] <- get_pubmed_distribution(node_DKIN)$Occurences
  wilkes_pcst_literatureBias_gamma[[as.character(g)]] <- get_pubmed_distribution(pcst_DKIN)$Occurences
}

wilkes_node_literatureBias_gamma_plot <- tibble(
  Values = c(wilkes_baselineKIN_df$Occurences, unname(unlist(wilkes_node_literatureBias_gamma))),
  Group = c(rep('Baseline KIN', nrow(wilkes_baselineKIN_df)), rep(names(wilkes_node_literatureBias_gamma), lengths(wilkes_node_literatureBias_gamma)))
)
wilkes_pcst_literatureBias_gamma_plot <- tibble(
  Values = c(wilkes_baselineKIN_df$Occurences, unname(unlist(wilkes_pcst_literatureBias_gamma))),
  Group = c(rep('Baseline KIN', nrow(wilkes_baselineKIN_df)), rep(names(wilkes_pcst_literatureBias_gamma), lengths(wilkes_pcst_literatureBias_gamma)))
)

# Convert Group labels to expressions for proper formatting
wilkes_node_literatureBias_gamma_plot$Group <- factor(
  wilkes_node_literatureBias_gamma_plot$Group,
  levels = c("Baseline KIN", names(wilkes_node_literatureBias_gamma)),
  labels = c("Baseline~KIN", paste0("gamma == ", names(wilkes_node_literatureBias_gamma)))
)
wilkes_pcst_literatureBias_gamma_plot$Group <- factor(
  wilkes_pcst_literatureBias_gamma_plot$Group,
  levels = c("Baseline KIN", names(wilkes_pcst_literatureBias_gamma)),
  labels = c("Baseline~KIN", paste0("gamma == ", names(wilkes_pcst_literatureBias_gamma)))
)

p1_gamma <- ggplot(wilkes_node_literatureBias_gamma_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Without PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  scale_fill_manual(
    values = c("steelblue", "tomato", "darkgreen", "purple", "orange")  # Adjust colors if needed
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

p2_gamma <- ggplot(wilkes_pcst_literatureBias_gamma_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "With PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  scale_fill_manual(
    values = c("steelblue", "tomato", "darkgreen", "purple", "orange")  # Adjust colors if needed
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

wilkes_node_literatureBias_beta <- list()
wilkes_pcst_literatureBias_beta <- list()
for (b in tested_beta) {
  node_DKIN <- read_tsv(paste0('../results/Hyperparameters/Wilkes2015/node_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma0_beta', b, '_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('../results/Hyperparameters/Wilkes2015/edge_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma0_beta', b, '_PCST_filtered_KIN.tsv'))
  
  wilkes_node_literatureBias_beta[[as.character(b)]] <- get_pubmed_distribution(node_DKIN)$Occurences
  wilkes_pcst_literatureBias_beta[[as.character(b)]] <- get_pubmed_distribution(pcst_DKIN)$Occurences
}

wilkes_node_literatureBias_beta_plot <- tibble(
  Values = c(wilkes_baselineKIN_df$Occurences, unname(unlist(wilkes_node_literatureBias_beta))),
  Group = c(rep('Baseline KIN', nrow(wilkes_baselineKIN_df)), rep(names(wilkes_node_literatureBias_beta), lengths(wilkes_node_literatureBias_beta)))
)
wilkes_pcst_literatureBias_beta_plot <- tibble(
  Values = c(wilkes_baselineKIN_df$Occurences, unname(unlist(wilkes_pcst_literatureBias_beta))),
  Group = c(rep('Baseline KIN', nrow(wilkes_baselineKIN_df)), rep(names(wilkes_pcst_literatureBias_beta), lengths(wilkes_pcst_literatureBias_beta)))
)

# Convert Group labels to expressions for proper formatting
wilkes_node_literatureBias_beta_plot$Group <- factor(
  wilkes_node_literatureBias_beta_plot$Group,
  levels = c("Baseline KIN", names(wilkes_node_literatureBias_beta)),
  labels = c("Baseline~KIN", paste0("beta == ", names(wilkes_node_literatureBias_beta)))
)
wilkes_pcst_literatureBias_beta_plot$Group <- factor(
  wilkes_pcst_literatureBias_beta_plot$Group,
  levels = c("Baseline KIN", names(wilkes_pcst_literatureBias_beta)),
  labels = c("Baseline~KIN", paste0("beta == ", names(wilkes_pcst_literatureBias_beta)))
)

p1_beta <- ggplot(wilkes_node_literatureBias_beta_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Without PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  #scale_fill_manual(
  #  values = c("steelblue", "tomato", "darkgreen", "purple")  # Adjust colors if needed
  #) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

p2_beta <- ggplot(wilkes_pcst_literatureBias_beta_plot, aes(x = Values, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "With PCST", x = "Number of PubMed mentions", y = "Density") +
  xlim(0, 1000) +
  #scale_fill_manual(
  #  values = c("steelblue", "tomato", "darkgreen", "purple")  # Adjust colors if needed
  #) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = scales::parse_format())

gamma_plot <- p1_gamma + 
  p2_gamma + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = expression(paste("Hyperparameter evaluation: Literature bias of ", gamma))) &
  theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm'))

beta_plot <- p1_beta + 
  p2_beta + 
  plot_layout(axis_titles = "collect", guides = "collect") +
  plot_annotation(title = expression(paste("Hyperparameter evaluation: Literature bias of ", beta))) &
  theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(0.5, 'cm'))
