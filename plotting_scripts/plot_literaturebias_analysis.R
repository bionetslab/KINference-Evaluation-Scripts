library(tidyverse)
library(ggplot2)
library(readxl)

biomart <- read_tsv('./data/biomart/mart_export_human_uniprot.txt')

gene2pubmed <- read_tsv('./data/gene2pubmed/gene2pubmed.tsv')
gene2ensembl <- read_tsv('./data/gene2pubmed/gene2ensembl.tsv')
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



bouhaddou_baselineKIN <- read_tsv('./results/Hyperparameters/Bouhaddou2023/baseline_KIN/Bouhaddou2023_n15_alpha0.9_KIN.tsv')
bouhaddou_baselineKIN_df <- get_pubmed_distribution(bouhaddou_baselineKIN)

wilkes_baseline_KIN <- read_tsv('./results/Hyperparameters/Wilkes2015/baseline_KIN/Wilkes2015_n15_alpha0.9_KIN.tsv')
wilkes_baselineKIN_df <- get_pubmed_distribution(wilkes_baseline_KIN)

tested_gamma <- c(0.5, 1, 1.5, 2, 2.5)
tested_beta <- c(0.2, 0.4, 0.6)
tested_delta <- c(0.7, 0.8, 0.9)


## BOUHADDOU et al.
bouhaddou_node_literatureBias_gamma <- list()
bouhaddou_pcst_literatureBias_gamma <- list()
for (g in tested_gamma) {
  node_DKIN <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/node_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma', g, '_beta0_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma', g, '_beta0_PCST_filtered_KIN.tsv'))
  
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

bouhaddou_gamma_data <- bind_rows(bouhaddou_node_literatureBias_gamma_plot %>% mutate(Type = 'Without PCST'), 
                                  bouhaddou_pcst_literatureBias_gamma_plot %>% mutate(Type = 'With PCST')) %>%
  mutate(Type = factor(Type, levels = c("Without PCST", "With PCST")))  # Ensure correct order

bouhaddou_gamma_plot <- ggplot(bouhaddou_gamma_data, aes(x = Group, y = Values, fill = Type)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  labs(title = expression(paste("DIFF")),
       x = "", y = "Number of PubMed mentions") +
  ylim(0, 500) +
  scale_x_discrete(labels = function(x) {parse(text = x)}) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    legend.title = element_blank(), 
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  guides(fill="none")

bouhaddou_node_literatureBias_beta <- list()
bouhaddou_pcst_literatureBias_beta <- list()
for (b in tested_beta) {
  node_DKIN <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/node_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma0_beta', b, '_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n15_alpha0.9_gamma0_beta', b, '_PCST_filtered_KIN.tsv'))
  
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

bouhaddou_beta_data <- bind_rows(bouhaddou_node_literatureBias_beta_plot %>% mutate(Type = 'Without PCST'), 
                                 bouhaddou_pcst_literatureBias_beta_plot %>% mutate(Type = 'With PCST')) %>%
  mutate(Type = factor(Type, levels = c("Without PCST", "With PCST")))  # Ensure correct order

bouhaddou_beta_plot <- ggplot(bouhaddou_beta_data, aes(x = Group, y = Values, fill = Type)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  labs(title = expression(paste("FS")),
       x = "", y = "Number of PubMed mentions") +
  ylim(0, 500) +
  scale_x_discrete(labels = function(x) {parse(text = x)}) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    legend.title = element_blank(), 
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  guides(fill="none")

bouhaddou_literatureBias_delta <- list()
for (d in tested_delta) {
  DKIN <- read_tsv(paste0('./results/Hyperparameters/Bouhaddou2023/edge_filtered_KIN/Bouhaddou2023_n10_alpha0.85_gamma0_beta0_delta', d, '_CORR_filtered_x0_KIN.tsv'))
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

bouhaddou_delta_data <- bouhaddou_literatureBias_delta_plot %>% mutate(Type = 'Without PCST') %>%
  mutate(Type = factor(Type, levels = c("Without PCST", "With PCST")))

bouhaddou_delta_plot <- ggplot(bouhaddou_delta_data, aes(x = Group, y = Values, fill = Type)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Dodge to separate conditions slightly
  labs(title = expression(paste("CORR")),
       x = "", y = "Number of PubMed mentions") +
  ylim(0, 500) +
  # Parse x-axis labels as expressions
  scale_x_discrete(labels = function(x) {parse(text = x)}) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    legend.title = element_blank(), 
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  guides(fill="none")

bouhaddou_plot <- ((bouhaddou_gamma_plot) | (bouhaddou_beta_plot) | (bouhaddou_delta_plot)) + 
  plot_layout(axis_titles = 'collect', guides = 'collect') +
  plot_annotation(tag_levels = list(c("A", "", "")), title = 'SARS-CoV-2') & 
  theme(plot.tag = element_text(face = "bold"), legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 18))


## WILKES et al.
wilkes_node_literatureBias_gamma <- list()
wilkes_pcst_literatureBias_gamma <- list()
for (g in tested_gamma) {
  node_DKIN <- read_tsv(paste0('./results/Hyperparameters/Wilkes2015/node_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma', g, '_beta0_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('./results/Hyperparameters/Wilkes2015/edge_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma', g, '_beta0_PCST_filtered_KIN.tsv'))
  
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

wilkes_gamma_data <- rbind(
  wilkes_node_literatureBias_gamma_plot %>% mutate(Type = 'Without PCST'), 
  wilkes_pcst_literatureBias_gamma_plot %>% mutate(Type = 'With PCST')
) %>% mutate(Type = factor(Type, levels = c('Without PCST', 'With PCST')))

wilkes_gamma_plot <- ggplot(wilkes_gamma_data, aes(x = Group, y = Values, fill = Type)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  labs(title = expression(paste("DIFF")),
       x = "", y = "Number of PubMed mentions") +
  ylim(0, 500) +
  scale_x_discrete(labels = function(x) {parse(text = x)}) +
  theme_minimal() +
  theme(
    legend.title = element_blank(), 
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

wilkes_node_literatureBias_beta <- list()
wilkes_pcst_literatureBias_beta <- list()
for (b in tested_beta) {
  node_DKIN <- read_tsv(paste0('./results/Hyperparameters/Wilkes2015/node_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma0_beta', b, '_KIN.tsv'))
  pcst_DKIN <- read_tsv(paste0('./results/Hyperparameters/Wilkes2015/edge_filtered_KIN/Wilkes2015_n15_alpha0.9_gamma0_beta', b, '_PCST_filtered_KIN.tsv'))
  
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

wilkes_beta_data <- rbind(
  wilkes_node_literatureBias_beta_plot %>% mutate(Type = 'Without PCST'),
  wilkes_pcst_literatureBias_beta_plot %>% mutate(Type = 'With PCST')
) %>% mutate(Type = factor(Type, levels = c('Without PCST', 'With PCST')))

wilkes_beta_plot <- ggplot(wilkes_beta_data, aes(x = Group, y = Values, fill = Type)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  labs(title = expression(paste("FS")),
       x = "", y = "Number of PubMed mentions") +
  ylim(0, 500) +
  scale_x_discrete(labels = function(x) {parse(text = x)}) +
  theme_minimal() +
  theme(
    legend.title = element_blank(), 
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

wilkes_plot <- (wilkes_gamma_plot + wilkes_beta_plot) + 
  plot_layout(axis_titles = 'collect', guides = 'collect') +
  plot_annotation(tag_levels = list("B"), title = 'PI3K inhibition in resistant breast cancer cell line') & 
  theme(plot.tag = element_text(face = "bold"), legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 18))



# perform correlation test:
print('------- SARS-CoV-2: -------')
bouhaddou_node_gamma_vals <- c(rep(0.5, length(bouhaddou_node_literatureBias_gamma$`0.5`)), rep(1.0, length(bouhaddou_node_literatureBias_gamma$`1`)), rep(1.5, length(bouhaddou_node_literatureBias_gamma$`1.5`)), rep(2.0, length(bouhaddou_node_literatureBias_gamma$`2`)), rep(2.5, length(bouhaddou_node_literatureBias_gamma$`2.5`)))
bouhaddou_pcst_gamma_vals <- c(rep(0.5, length(bouhaddou_pcst_literatureBias_gamma$`0.5`)), rep(1.0, length(bouhaddou_pcst_literatureBias_gamma$`1`)), rep(1.5, length(bouhaddou_pcst_literatureBias_gamma$`1.5`)), rep(2.0, length(bouhaddou_pcst_literatureBias_gamma$`2`)), rep(2.5, length(bouhaddou_pcst_literatureBias_gamma$`2.5`)))
bouhaddou_node_beta_vals <- c(rep(1, length(bouhaddou_node_literatureBias_beta$`0.2`)), rep(2, length(bouhaddou_node_literatureBias_beta$`0.4`)), rep(3, length(bouhaddou_node_literatureBias_beta$`0.6`)))
bouhaddou_pcst_beta_vals <- c(rep(1, length(bouhaddou_pcst_literatureBias_beta$`0.2`)), rep(2, length(bouhaddou_pcst_literatureBias_beta$`0.4`)), rep(3, length(bouhaddou_pcst_literatureBias_beta$`0.6`)))
bouhaddou_node_delta_vals <- c(rep(0.7, length(bouhaddou_literatureBias_delta$`0.7`)), rep(0.8, length(bouhaddou_literatureBias_delta$`0.8`)), rep(0.9, length(bouhaddou_literatureBias_delta$`0.9`)))
print('------- Without PCST: -------')
print('------- Gamma: -------')
cor.test(bouhaddou_node_gamma_vals, unname(unlist(bouhaddou_node_literatureBias_gamma)))
print('------- Beta: -------')
cor.test(bouhaddou_node_beta_vals, unname(unlist(bouhaddou_node_literatureBias_beta)))
print('------- Delta: -------')
cor.test(bouhaddou_node_delta_vals, unname(unlist(bouhaddou_literatureBias_delta)))

print('With PCST:')
print('------- Gamma: -------')
cor.test(bouhaddou_pcst_gamma_vals, unname(unlist(bouhaddou_pcst_literatureBias_gamma)))
print('------- Beta: -------')
cor.test(bouhaddou_pcst_beta_vals, unname(unlist(bouhaddou_pcst_literatureBias_beta)))


print('------- PI3K Inhibitor: -------')
wilkes_node_gamma_vals <- c(rep(0.5, length(wilkes_node_literatureBias_gamma$`0.5`)), rep(1.0, length(wilkes_node_literatureBias_gamma$`1`)), rep(1.5, length(wilkes_node_literatureBias_gamma$`1.5`)), rep(2.0, length(wilkes_node_literatureBias_gamma$`2`)), rep(2.5, length(wilkes_node_literatureBias_gamma$`2.5`)))
wilkes_pcst_gamma_vals <- c(rep(0.5, length(wilkes_pcst_literatureBias_gamma$`0.5`)), rep(1.0, length(wilkes_pcst_literatureBias_gamma$`1`)), rep(1.5, length(wilkes_pcst_literatureBias_gamma$`1.5`)), rep(2.0, length(wilkes_pcst_literatureBias_gamma$`2`)), rep(2.5, length(wilkes_pcst_literatureBias_gamma$`2.5`)))
wilkes_node_beta_vals <- c(rep(1, length(wilkes_node_literatureBias_beta$`0.2`)), rep(2, length(wilkes_node_literatureBias_beta$`0.4`)), rep(3, length(wilkes_node_literatureBias_beta$`0.6`)))
wilkes_pcst_beta_vals <- c(rep(1, length(wilkes_pcst_literatureBias_beta$`0.2`)), rep(2, length(wilkes_pcst_literatureBias_beta$`0.4`)), rep(3, length(wilkes_pcst_literatureBias_beta$`0.6`)))
print('------- Without PCST: -------')
print('------- Gamma: -------')
cor.test(wilkes_node_gamma_vals, unname(unlist(wilkes_node_literatureBias_gamma)))
print('------- Beta: -------')
cor.test(wilkes_node_beta_vals, unname(unlist(wilkes_node_literatureBias_beta)))

print('With PCST:')
print('------- Gamma: -------')
cor.test(wilkes_pcst_gamma_vals, unname(unlist(wilkes_pcst_literatureBias_gamma)))
print('------- Beta: -------')
cor.test(wilkes_pcst_beta_vals, unname(unlist(wilkes_pcst_literatureBias_beta)))



#####################
## Invergo et al. ##
####################
invergo_baseline_KIN <- read_tsv('./results/baselineKIN_comparisons/baseline_KIN/invergo2020_comparison_KIN.tsv')
invergo_baseline_KIN_interactions <- unique(paste0(invergo_baseline_KIN$Source, '_', invergo_baseline_KIN$Target_Uniprot))

invergo_baseline_KIN_without_n_threshold <- read_tsv('./results/baselineKIN_comparisons_without_n_threshold/baseline_KIN/invergo2020_comparison_KIN.tsv') 

# get invergo interactions and map them to UniProt IDs
invergo_data <- read_excel('./data/competitors/invergo2020_interactions.xlsx', na = c('NA', '')) 
colnames(invergo_data) <- invergo_data[1, ]
invergo_data <- invergo_data[-1, ]
# filtering based on threshold in invergo paper
invergo_data <- invergo_data %>% mutate(Mean_posterior_probability_of_regulation = as.double(Mean_posterior_probability_of_regulation)) %>% filter(Mean_posterior_probability_of_regulation > 0.5, !is.na(Mean_posterior_probability_of_regulation))
invergo_data$Source <- sapply(
  invergo_data$Transducing_kinase, 
  function(x) {unname(uniprots[x])[1]}
)
invergo_data$Target <- sapply(
  invergo_data$Substrate_kinase, 
  function(x) {unname(uniprots[x])[1]}
)
get_pubmed_distribution_invergo <- function(data) {
  targets <- unique(data$Substrate_kinase)
  targets <- targets[!is.na(targets)]
  ensg_IDs <- unique(biomart[which(biomart$`UniProtKB Gene Name symbol` %in% targets), ][['Gene stable ID']])
  geneIDs <- unique(gene2ensembl[which(gene2ensembl$Ensembl_gene_identifier %in% ensg_IDs),]$'GeneID')
  distribution <- table(gene2pubmed[which(gene2pubmed$GeneID %in% geneIDs),]$GeneID)
  df <- tibble(
    "Gene" = names(distribution),
    "Occurences" = as.integer(unname(distribution))
  )
  return(df)
}
invergo_lit_bias <- get_pubmed_distribution_invergo(invergo_data)

comp_data <- tibble(
  Values = invergo_lit_bias$Occurences,
  Group = rep("Invergo et al.", nrow(invergo_lit_bias))
)

###################
## Buljan et al. ##
###################
buljan_baseline_KIN <- read_tsv('./results/baselineKIN_comparisons/baseline_KIN/buljan2020_comparison_KIN.tsv')

buljan_data <- read_excel('./data/competitors/buljan2020_interactions.xlsx', na = c('NA', ''))
colnames(buljan_data) <- buljan_data[2,]
buljan_data <- buljan_data[-c(1, 2),]
# filtering based on threshold in buljan paper
buljan_data <- buljan_data %>% mutate(wdscore = as.double(wdscore), gfpratio = as.double(gfpratio)) %>% filter(wdscore > 73.6, gfpratio > 18.4)
get_pubmed_distribution_buljan <- function(data) {
  targets <- unique(data$Protein_id)
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
buljan_lit_bias <- get_pubmed_distribution_buljan(buljan_data) 

comp_data <- rbind(comp_data,
  tibble(
    Values = buljan_lit_bias$Occurences,
    Group = rep("Buljan et al.", nrow(buljan_lit_bias))
  )
)

###################
## GRNBoost2 ##
###################
grnboost2_x0 <- read_tsv('./results/GRNBoost2/x0_network.tsv', col_names = F)
grnboost2_x1 <- read_tsv('./results/GRNBoost2/x1_network.tsv', col_names = F)
grnboost2_x0$Source <- sapply(strsplit(grnboost2_x0$X1, '_'), function(x){ x[1] })
grnboost2_x1$Source <- sapply(strsplit(grnboost2_x1$X1, '_'), function(x){ x[1] })
grnboost2_x0$Target <- sapply(strsplit(grnboost2_x0$X2, '_'), function(x){ x[1] })
grnboost2_x1$Target <- sapply(strsplit(grnboost2_x1$X2, '_'), function(x){ x[1] })
grnboost2_data <- tibble(Source = c(grnboost2_x0$Source, grnboost2_x1$Source), Target = c(grnboost2_x0$Target, grnboost2_x1$Target))
grnboost2_data <- unique(grnboost2_data)

get_pubmed_distribution_grnboost2 <- function(data) {
  targets <- unique(data$Target)
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

grnboost2_lit_bias <- get_pubmed_distribution_grnboost2(grnboost2_data) 

comp_data <- rbind(comp_data,
                   tibble(
                     Values = grnboost2_lit_bias$Occurences,
                     Group = rep("GRNBoost2", nrow(grnboost2_lit_bias))
                   )
)

###################
## KINference ##
###################
baseline_KIN <- read_tsv('./results/example_results_baseline_KIN/example_matrix_run_KIN.tsv')
baseline_KIN <- baseline_KIN %>% filter(Edge_type == 'KS')
baseline_KIN_lit_bias <- get_pubmed_distribution(baseline_KIN)
comp_data <- rbind(comp_data,
                   tibble(
                     Values = baseline_KIN_lit_bias$Occurences,
                     Group = rep("KINference", nrow(baseline_KIN_lit_bias))
                   )
)

# PLOT
comp_data$Group <- factor(
  comp_data$Group,
  levels = c("Invergo et al.", "Buljan et al.", "GRNBoost2", "KINference"),
  labels = c("Invergo~et~al.", "Buljan~et~al.", "GRNBoost2", "KINference")
)
competitor_plot <- ggplot(comp_data, aes(x = Group, y = Values, fill = Group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  labs(title = expression(paste("Literature bias: KINference vs. Competitors")),
       x = "", y = "Number of PubMed mentions") +
  ylim(0, 500) +
  scale_x_discrete(labels = function(x) {parse(text = x)}) +
  theme_minimal() +
  theme(
    legend.title = element_blank(), 
    axis.text = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none"
  )