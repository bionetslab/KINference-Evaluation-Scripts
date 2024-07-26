library(tidyverse)
library(data.table)
library(UniProt.ws)
library(funscoR)
library(OmnipathR)
library(argparse)

### Helper functions
# Uniprot sequence loading helper function
.query_uniprot_ws_for_sequence <- function(uniprots, fields = c("accession","sequence")){
  queries <- paste0("accession:", uniprots)
  res <- setDT(UniProt.ws::queryUniProt(queries, fields = fields))
  #print (res)
  res[uniprots,query := Entry , on = "Entry"]
  res[]
  
}
# Uniprot sequence loading helper function
uniprotSequencesFromWeb <- function (uniprotIDs, chunkSize = 25, fields = c("accession","sequence")){
  uniqueUniprots <- unique(uniprotIDs)
  chunks <- split(uniqueUniprots, (1:length(uniqueUniprots))%/%chunkSize)
  
  infoMapList <- pbapply::pblapply(chunks, .query_uniprot_ws_for_sequence, fields = fields)
  
  infoMap <- rbindlist(infoMapList)
  return (infoMap)
}


#' Load kinase data for substrate scoring
#' 
#' This function loads the kinase motifs, apriori distributions, scaling factors, and name mappings
#' required for substrate scoring.
#' 
#' @param kinase_motifs.path Path to the file containing kinase motifs data.
#' @param kinase_aprior_distributions.path Path to the directory containing apriori distributions for kinases.
#' @param kinase_scaling_factors.path Path to the file containing scaling factors for substrate scoring.
#' @param kinase_name_mappings.path Path to the file containing kinase name mappings for human.
#' @param kinase_name_mappings_mouse.path Path to the file containing kinase name mappings for mouse.
#' @return A list containing the loaded data.
#' @examples
#' # Example usage
#' load_kinase_data()
load_kinase_data <- function(
  kinase_motifs.path,
  kinase_aprior_distributions.path,
  kinase_name_mappings.path) {
  
  # reading in kinase motif data
  kinase_motifs <- fread(kinase_motifs.path)

  # reading in mapping of kinases to uniprot id file
  kinase_names_to_uniprotmapping <- fread(kinase_name_mappings.path)

  kinases <- kinase_names_to_uniprotmapping[['ACC#']]

  # Computing a-priori distributions of kinases based on log2 score database
  kinase_aprior_distributions <- list()
  for (i in seq_len(nrow(kinase_motifs))) {
    row <- kinase_motifs[i, ]
    kin <- row$KINASE
    kin_aprior_file_name <- paste(kinase_aprior_distributions.path, kinase_motifs$KINASE[i], '.txt', sep = '')
    aprior_distribution <- ecdf(fread(kin_aprior_file_name)[[1]])
    kinase_aprior_distributions <- append(kinase_aprior_distributions, aprior_distribution)
  }
  names(kinase_aprior_distributions) <- kinases
  
  return(list(
    kinase_motifs = kinase_motifs, 
    kinase_aprior_distributions = kinase_aprior_distributions, 
    kinase_name_mappings = kinase_names_to_uniprotmapping
  ))

}


#' Perform substrate scoring
#' 
#' This function calculates substrate scores as described in https://doi.org/10.1038/s41586-022-05575-3 Fig. 3
#' 
#' @param intensities A data frame containing intensities.
#' @param output.path Path to the directory where the output file will be saved.
#' @param output.id Name of the output file.
#' @param kinase_data Kinase information loaded by the function load_kinase_data()
#' @param percentile_rank_threshold Threshold for percentile rank (default: 15).
#' @param alpha Threshold for percentile scores based of the percentage (default: 0.9).
#' @return A data frame containing the calculated substrate scores.
#' @examples
#' # Example usage
#' substrate_scoring(intensities_df, "./output/", "cond1_vs_cond2.tsv", kinase_motifs_df, kinase_apriori_distributions_list, kinase_scaling_factors_df, 15)
infer_baseline_tyrosine_KIN <- function(
    md, output.path, output.id, kinase_data, 
    n = 15, alpha = 0.9) {

  kinases <- kinase_data$kinase_name_mappings$`ACC#`
  # Loading sequences from Uniprot  
  sequences <- uniprotSequencesFromWeb(unique(md$Uniprot))

  get_modified_sequence <- function(seq, pos) {
    if (pos > nchar(seq) || !(str_sub(seq, pos, pos) %in% c('Y'))) {
      return(NULL)
    } else if ((pos > 5) && ((pos + 5) < nchar(seq))) {
      return(list(seq = str_sub(seq, pos - 5, pos + 5), zero_idx = 6))
    } else if (pos > 5) {
      return(list(seq = str_sub(seq, pos - 5), zero_idx = 6))
    } else if ((pos + 5) < nchar(seq)) {
      n_seq <- str_sub(seq, 1, pos + 5)
      return(list(seq = n_seq, zero_idx = nchar(n_seq) - 4))
    } else {
      return(NULL)
    }
  }
  
  sequences_list <- lapply(1:nrow(md), function(i) {
    seq <- sequences$Sequence[which(sequences$query == md$Uniprot[i])]
    pos <- md$Pos[i]
    get_modified_sequence(seq, pos)
  })
  
  valid_indices <- which(!sapply(sequences_list, is.null))
  sequences_list <- sequences_list[valid_indices]
  md_filtered <- md[valid_indices,]
  
  get_scores_log2 <- function(seq_info, kinase_data) {
    seq <- seq_info$seq
    zero_idx <- seq_info$zero_idx
    seq_split <- str_extract_all(seq, stringr::boundary('character'))[[1]]
    indices <- if (zero_idx == 6) {
      seq(-5, length(seq_split) - 6, 1)
    } else {
      seq(-(length(seq_split) - 5), 5, 1)
    }
    
    aa <- paste0(indices, seq_split)[-zero_idx]
    scores_log2 <- log2(apply(kinase_data$kinase_motifs[, ..aa], 1, prod))
    names(scores_log2) <- kinases
    scores_log2
  }
  
  scores_log2_list <- lapply(sequences_list, get_scores_log2, kinase_data = kinase_data)
  
  percentile_scores_list <- lapply(scores_log2_list, function(scores_log2) {
    sapply(names(kinase_data$kinase_aprior_distributions), function(x) { 
      kinase_data$kinase_aprior_distributions[[x]](scores_log2[[x]]) 
    })
  })
  
  construct_final_dataframe <- function(i, scores_log2, percentile_scores) {
    percentile_scores <- percentile_scores[order(-percentile_scores)]
    edges_threshold <- min(n, length(which(percentile_scores >= alpha)))
    if (edges_threshold == 0) {
      return(NULL)
    }
    
    seq_info <- sequences_list[[i]]
    new_rows <- data.frame(
      Source = names(percentile_scores)[1:edges_threshold], 
      Target = rep(md_filtered$Protein[i], edges_threshold), 
      Target_Uniprot = rep(md_filtered$Uniprot[i], edges_threshold),
      ModifiedSequence = rep(seq_info$seq, edges_threshold), 
      log2Score = scores_log2[names(percentile_scores)[1:edges_threshold]], 
      percentileScore = percentile_scores[1:edges_threshold], 
      percentileRank = seq(1, edges_threshold, 1), 
      f = rep(md_filtered$f[i], edges_threshold)
    )
    return(new_rows)
  }
  
  result_list <- mapply(construct_final_dataframe, 
                        i = seq_along(scores_log2_list), 
                        scores_log2 = scores_log2_list, 
                        percentile_scores = percentile_scores_list, 
                        SIMPLIFY = FALSE)
  
  result_list <- result_list[which(!sapply(result_list, is.null))]
  baseline_KIN <- rbindlist(result_list)
 
  baseline_KIN_KS_edges <- baseline_KIN
  baseline_KIN_KS_edges$Type = 'KS'
  SK_edges <- tibble(
    Source =  unique(baseline_KIN_KS_edges$Target[which(baseline_KIN_KS_edges$Target_Uniprot %in% kinase_data$kinase_name_mappings$`ACC#`)]),
    Target = sapply(strsplit(Source, '_'), function(x) { x[1] }),
    Target_Uniprot = Target,
    ModifiedSequence = '',
    log2Score = 0,
    percentileScore = 0,
    percentileRank = 0,
    f = 0,
    Type = 'SK'
  )
  
  write_tsv(
    rbind(baseline_KIN_KS_edges, SK_edges),
    paste0(output.path, '/baseline_KIN/', output.id, '_tyrosine_kinases.tsv')
  )
  
  return(baseline_KIN_KS_edges)
}

infer_baseline_serine_threonine_KIN <- function(
    md, output.path, output.id, kinase_data, 
    n = 15, alpha = 0.9) {

  kinases <- kinase_data$kinase_name_mappings$`ACC#`
  # Loading sequences from Uniprot  
  sequences <- uniprotSequencesFromWeb(unique(md$Uniprot))

  # Create output directory
  dir.create(paste0(output.path, '/baseline_KIN/'), recursive = TRUE)

  get_modified_sequence <- function(seq, pos) {
    if (pos > nchar(seq) || !(str_sub(seq, pos, pos) %in% c('S', 'T'))) {
      return(NULL)
    } else if ((pos > 5) && ((pos + 4) < nchar(seq))) {
      return(list(seq = str_sub(seq, pos - 5, pos + 4), zero_idx = 6))
    } else if (pos > 5) {
      return(list(seq = str_sub(seq, pos - 5), zero_idx = 6))
    } else if ((pos + 4) < nchar(seq)) {
      n_seq <- str_sub(seq, 1, pos + 4)
      return(list(seq = n_seq, zero_idx = nchar(n_seq) - 4))
    } else {
      return(NULL)
    }
  }
  
  sequences_list <- lapply(1:nrow(md), function(i) {
    seq <- sequences$Sequence[which(sequences$query == md$Uniprot[i])]
    pos <- md$Pos[i]
    get_modified_sequence(seq, pos)
  })
  
  valid_indices <- which(!sapply(sequences_list, is.null))
  sequences_list <- sequences_list[valid_indices]
  md_filtered <- md[valid_indices,]
  
  get_scores_log2 <- function(seq_info, kinase_data) {
    seq <- seq_info$seq
    zero_idx <- seq_info$zero_idx
    seq_split <- str_extract_all(seq, stringr::boundary('character'))[[1]]
    indices <- if (zero_idx == 6) {
      seq(-5, length(seq_split) - 6, 1)
    } else {
      seq(-(length(seq_split) - 5), 4, 1)
    }
    
    aa <- paste0(indices, seq_split)[-zero_idx]
    scores_log2 <- log2(apply(kinase_data$kinase_motifs[, ..aa], 1, prod))
    names(scores_log2) <- kinases
    scores_log2
  }
  
  scores_log2_list <- lapply(sequences_list, get_scores_log2, kinase_data = kinase_data)
  
  percentile_scores_list <- lapply(scores_log2_list, function(scores_log2) {
    sapply(names(kinase_data$kinase_aprior_distributions), function(x) { 
      kinase_data$kinase_aprior_distributions[[x]](scores_log2[[x]]) 
    })
  })
  
  construct_final_dataframe <- function(i, scores_log2, percentile_scores) {
    percentile_scores <- percentile_scores[order(-percentile_scores)]
    edges_threshold <- min(n, length(which(percentile_scores >= alpha)))
    if (edges_threshold == 0) {
      return(NULL)
    }
    
    seq_info <- sequences_list[[i]]
    new_rows <- data.frame(
      Source = names(percentile_scores)[1:edges_threshold], 
      Target = rep(md_filtered$Protein[i], edges_threshold), 
      Target_Uniprot = rep(md_filtered$Uniprot[i], edges_threshold),
      ModifiedSequence = rep(seq_info$seq, edges_threshold), 
      log2Score = scores_log2[names(percentile_scores)[1:edges_threshold]], 
      percentileScore = percentile_scores[1:edges_threshold], 
      percentileRank = seq(1, edges_threshold, 1), 
      f = rep(md_filtered$f[i], edges_threshold)
    )
    return(new_rows)
  }
  
  result_list <- mapply(construct_final_dataframe, 
                        i = seq_along(scores_log2_list), 
                        scores_log2 = scores_log2_list, 
                        percentile_scores = percentile_scores_list, 
                        SIMPLIFY = FALSE)
  
  result_list <- result_list[which(!sapply(result_list, is.null))]
  baseline_KIN <- rbindlist(result_list)
 
  baseline_KIN_KS_edges <- baseline_KIN
  baseline_KIN_KS_edges$Type = 'KS'
  SK_edges <- tibble(
    Source =  unique(baseline_KIN_KS_edges$Target[which(baseline_KIN_KS_edges$Target_Uniprot %in% kinase_data$kinase_name_mappings$`ACC#`)]),
    Target = sapply(strsplit(Source, '_'), function(x) { x[1] }),
    Target_Uniprot = Target,
    ModifiedSequence = '',
    log2Score = 0,
    percentileScore = 0,
    percentileRank = 0,
    f = 0,
    Type = 'SK'
  )
  
  write_tsv(
    rbind(baseline_KIN_KS_edges, SK_edges),
    paste0(output.path, '/baseline_KIN/', output.id, '_serine_threonine_kinases.tsv')
  )
  
  return(baseline_KIN_KS_edges)
}


#' Perform kinase enrichment analysis
#' 
#' This function performs kinase enrichment analysis based on intensities and substrate scores,
#' and writes the results to a specified output file.
#' 
#' @param intensities A data frame containing intensities.
#' @param sub_scores A data frame containing substrate scores.
#' @param output.path Path to the directory where the output file will be saved.
#' @param output.id Name of the output file.
#' @param gamma Threshold for log2 fold change.
#' @param kinase_data Kinase information loaded by the function load_kinase_data()
#' @return A data frame containing the results of the kinase enrichment analysis.
#' @examples
#' # Example usage
#' kinase_enrichment(intensities_df, substrate_scores_df, "./output/", "cond1_vs_cond2.tsv", 1.0)
kinase_enrichment <- function(
    md, baseline_KIN, output.path, output.id, kinase_data, gamma = 1.0) {
  
  kinases <- kinase_data$kinase_name_mappings$`ACC#`
  # creating output directory for enrichments
  dir.create(paste0(output.path, '/enrichments/'), recursive = TRUE)
  
  # Upregulated set := Phosphosites with log2FC > gamma
  # Downregulated set := Phosphosites with log2FC < - (gamma)
  # Background set := Phosphosites with -(gamma) <= log2FC <= gamma
  upregulated_phosphorylationSites <- md[which(md$f > gamma),]
  downregulated_phosphorylationSites <- md[which(md$f < -gamma),]
  background_phosphorylationSites <- md[which(abs(md$f) <= gamma),]

  # Getting all interactions in each set 
  upregulated_interactions <- baseline_KIN[which(baseline_KIN$Target %in% upregulated_phosphorylationSites$Protein),]
  downregulated_interactions <- baseline_KIN[which(baseline_KIN$Target %in% downregulated_phosphorylationSites$Protein),]
  background_interactions <- baseline_KIN[which(baseline_KIN$Target %in% background_phosphorylationSites$Protein),]

  # Counting the occurrences of each kinase in the three sets with 
  up_counts <- data.frame(table(upregulated_interactions$Source))
  if (nrow(up_counts) != length(kinases)) up_counts <- rbind(up_counts, data.frame(Var1=setdiff(kinases, up_counts$Var1), Freq=0))
  down_counts <- data.frame(table(downregulated_interactions$Source))
  if (nrow(down_counts) != length(kinases)) down_counts <- rbind(down_counts, data.frame(Var1=setdiff(kinases, down_counts$Var1), Freq=0))
  background_counts <- data.frame(table(background_interactions$Source))
  if (nrow(background_counts) != length(kinases)) background_counts <- rbind(background_counts, data.frame(Var1=setdiff(kinases, background_counts$Var1), Freq=0))

  # Enrichment scoring
  enrich_val_up <- c()
  enrich_val_down <- c()
  dominant_enrichment_value <- c()
  dominant_enrichment_log2_value <- c()
  p_val_up <- c()
  p_val_down <- c()
  dominant_p_value <- c()
  dominant_directions <- c()

  up_set_size <- nrow(upregulated_interactions)
  down_set_size <- nrow(downregulated_interactions)
  background_set_size <- nrow(background_interactions)
  for (i in 1:nrow(up_counts)) {

    up_hits <- up_counts$Freq[i]
    down_hits <- down_counts$Freq[i]
    background_hits <- background_counts$Freq[i]
    
    #TODO: Fix Haldane correction? Adding 0.5 to every cell of the contigency table in R does not seem to work for fisher.test
    # Haldane correction up_hits (No kinase hits in the upregulated (or background) interactions)
    if (up_hits == 0 | background_hits == 0) {
      d_up <- data.frame(
        x = c(up_hits + 0.5, background_hits + 0.5),
        y = c(up_set_size + 0.5, background_set_size + 0.5)
      )
    } else {
      d_up <- data.frame(
        x = c(up_hits, background_hits),
        y = c(up_set_size, background_set_size)
      )
    }
    
    # Haldane correction down_hits
    if (down_hits == 0 | background_hits == 0) {
      d_down <- data.frame(
        x = c(down_hits + 0.5, background_hits + 0.5),
        y = c(down_set_size + 0.5, background_set_size + 0.5)
      )
    } else {
      d_down <- data.frame(
        x = c(down_hits, background_hits),
        y = c(down_set_size, background_set_size)
      )
    }
    
    # Calculating p-value and enrichment score ((up_hits/up_set_size)/(background_hits/background_set_size))
    # using one sided fisher exact z test (assumption: background set is bigger than up-,downregulatory set)
    # Dominant value: max(enrich_up, enrich_down)
    test_up <- fisher.test(d_up, alternative = 'greater')
    test_down <- fisher.test(d_down, alternative = 'greater')
    
    enrich_up <- unname(test_up[['estimate']])
    enrich_down <- unname(test_down[['estimate']])
    
    enrich_val_up <- c(enrich_val_up, enrich_up)
    enrich_val_down <- c(enrich_val_down, enrich_down)
    dominant_value <- ifelse(enrich_up > enrich_down, 1, 2)
    dominant_directions <- c(dominant_directions, c('upregulated set', 'downregulated set')[dominant_value])
    dominant_enrichment_log2_value <- c(dominant_enrichment_log2_value, ifelse(enrich_up > enrich_down, log2(enrich_up), -log2(enrich_down)))
    dominant_enrichment_value <- c(dominant_enrichment_value, ifelse(enrich_up > enrich_down, enrich_up, 2 ** -log2(enrich_down)))
    p_val_up <- c(p_val_up, test_up[['p.value']])
    p_val_down <- c(p_val_down, test_down[['p.value']])
    dominant_p_value <- c(dominant_p_value, c(test_up[['p.value']], test_down[['p.value']])[dominant_value])
  }

  # Adjusting p values using Benjamini Hochberg correction
  p_val_up_adjusted <- p.adjust(p_val_up, method = 'BH')
  p_val_down_adjusted <- p.adjust(p_val_down, method = 'BH')
  dominant_p_value_adjusted <- p.adjust(dominant_p_value, method = 'BH')

  result <- data.frame(uniprotID = up_counts$Var1, 
                      gene = sapply(up_counts$Var1, function(x) {kinase_data$kinase_name_mappings$GENE[which(kinase_data$kinase_name_mappings$`ACC#` == x)]}),
                      upregulated_set_hits = up_counts$Freq, upregulated_set_size = up_set_size, 
                      downregulated_set_hits = down_counts$Freq, down_regulated_set_size = down_set_size, 
                      background_set_hits = background_counts$Freq, background_set_size = background_set_size, 
                      upregulated_enrichment_value = enrich_val_up, upregulated_enrichment_value_log2 = log2(enrich_val_up), 
                      upregulated_p_value = p_val_up, upregulated_adjusted_p_value = p_val_up_adjusted, upregulated_p_value_log10_abs = abs(log10(p_val_up)), 
                      upregulated_adjusted_p_value_log10_abs = abs(log10(p_val_up_adjusted)), downregulated_enrichment_value = enrich_val_down,
                      downregulated_enrichment_value_log2 = log2(enrich_val_down), downregulated_p_value = p_val_down,
                      downregulated_adjusted_p_value = p_val_down_adjusted, downregulated_p_value_log10_abs = abs(log10(p_val_down)), downregulated_adjusted_p_value_log10_abs = abs(log10(p_val_down_adjusted)),
                      dominant_enrichment_value = dominant_enrichment_value, dominant_enrichment_value_log2 = dominant_enrichment_log2_value,
                      dominant_p_value = dominant_p_value, dominant_adjusted_p_value = dominant_p_value_adjusted,
                      dominant_p_value_log10_abs = abs(log10(dominant_p_value)), dominant_adjusted_p_value_log10_abs = abs(log10(dominant_p_value_adjusted)),
                      dominant_direction = dominant_directions
                    )

  write_tsv(result, paste0(output.path, '/enrichments/', output.id, '.tsv'))
}


#' Perform functional scoring (!!! Only works for human data)
#' 
#' This function performs functional scoring based on intensities and writes the results
#' to a specified output file (see https://doi.org/10.1038/s41587- 019-0344-3).
#' 
#' @param intensities A data frame containing intensities.
#' @param output.path Path to the directory where the output file will be saved.
#' @return A data frame containing the results of the functional scoring.
#' @examples
#' # Example usage
#' functional_scoring(intensities_df, "./output/")
compute_FS_filter <- function(md, output.path, output.id) {
  
  data_targets <- data.frame(
    acc=md$Uniprot,
    position=md$Pos,
    residue=md$AA
  )
  
  merged_phosphoproteome <- unique(rbind(data_targets, phosphoproteome))
  
  # removing duplicates (keeping reads from lab data)
  dup <- duplicated(merged_phosphoproteome[, c(1,2)])
  merged_phosphoproteome <- merged_phosphoproteome[!dup,]
  
  annotated_phos <- annotate_sites(merged_phosphoproteome)
  
  ST_features <- preprocess_features(annotated_phos, "ST")
  Y_features <- preprocess_features(annotated_phos, "Y")
  
  ## train new model
  ST_model <- train_funscore(ST_features, "ST", psp, ncores = 8)
  Y_model <- train_funscore(Y_features, "Y", psp, ncores = 8)
  
  ## predict funcscoR for all sites
  ST_scores <- predict_funscore(ST_features, ST_model, ncores = 8)
  Y_scores <- predict_funscore(Y_features, Y_model, ncores = 8)
  
  ## gather all predictions
  all_scores <- bind_rows(ST_scores, Y_scores) %>%
    mutate(probabilities = log_scaling(probabilities))
  
  write_tsv(all_scores, paste0(output.path, '/', output.id, '_FSscores.tsv'))
  
  return(all_scores)
}


#' Perform correlation scoring
#' 
#' This function calculates the correlation between intensities and substrate scores
#' and writes the results to a specified output file.
#' 
#' @param intensities A data frame containing intensities.
#' @param sub_scores A data frame containing substrate scores.
#' @param condition Condition in intensity dataframe that is being used for correlation scoring
#' @param output.path Path to the directory where the output file will be saved.
#' @param output.id Name of the output file.
#' @param output.condition Condition that is being currently scored
#' @param m Threshold for the number of samples needed for correlation calculation (default: 10).
#' @param delta Threshold for correlation coefficient
#' @param epsilon Threshold for significance of correlation test
#' @return Edge list E_CORR
compute_CORR_filter <- function(
  x, md, baseline_KIN, output.path, output.id, serine_threonine_kinase_data,
  tyrosine_kinase_data, m = 10, delta = 0.8, epsilon = 0.05) {
  
  dir.create(paste0(output.path, '/correlation_networks/'), recursive = TRUE)
  
  indices <- unname(which((ncol(x) - rowSums(is.na(x))) >= m))
  x <- x[indices, ]
  md <- md[indices, ]
  
  # filter out all downstream interactions that have lower than m intensity measurements
  baseline_KIN <- baseline_KIN[which(baseline_KIN$Target %in% md$Protein),]

  # only keeping the edges that have an intensity measure on the source (the inferred network from the substrate scoring is kinase -> PhosphoSite not PhosphoSite -> PhosphoSite)
  baseline_KIN <- baseline_KIN[baseline_KIN$Source %in% md$Uniprot, ]
  
  keep_edge <- rep(F, nrow(baseline_KIN))
  
  for (kinase in unique(baseline_KIN$Source)) {
    edge_indices <- which(baseline_KIN$Source == kinase)
    kinase_phosphorylationSites <- unique(baseline_KIN$Target[which(baseline_KIN$Target_Uniprot == kinase)])
    for (index in edge_indices) {
      target <- baseline_KIN$Target[index]
      for (source in kinase_phosphorylationSites) {
        corr_test <- cor.test(unlist(x[which(rownames(x) == source),]), unlist(x[which(rownames(x) == target),]))
        if (corr_test$p.value <= epsilon & abs(corr_test$estimate) >= delta) {
          keep_edge[index] <- T
          break
        }
      }
    }
  }
  
  corr_KIN <- baseline_KIN[keep_edge,]
  non_targeted_kinases <- unique(baseline_KIN$Source[which(!(baseline_KIN$Source %in% baseline_KIN$Target_Uniprot))])
  
  corr_KIN <- rbind(corr_KIN, baseline_KIN[which(baseline_KIN$Source %in% non_targeted_kinases),])
  
  corr_KIN_KS_edges <- baseline_KIN
  corr_KIN_KS_edges$Type = 'KS'
  SK_edges <- tibble(
    Source =  unique(corr_KIN_KS_edges$Target[which(corr_KIN_KS_edges$Target_Uniprot %in% c(serine_threonine_kinase_data$kinase_name_mappings$`ACC#`, tyrosine_kinase_data$kinase_name_mappings$`ACC#`))]),
    Target = sapply(strsplit(Source, '_'), function(x) { x[1] }),
    Target_Uniprot = Target,
    ModifiedSequence = '',
    log2Score = 0,
    percentileScore = 0,
    percentileRank = 0,
    f = 0,
    FS = 0,
    Type = 'SK'
  )
  
  write_tsv(
    rbind(corr_KIN_KS_edges, SK_edges),
    paste0(output.path, '/correlation_networks/', output.id, '.tsv')
  )
  #write_tsv(corr_KIN, paste0(output.path, '/correlation_networks/', output.id, '.tsv'))
  
  return(rbind(corr_KIN_KS_edges, SK_edges))
}

#' Run the PCST algorithm on the given baseline KIN
#' 
#' @param sub_scores Baseline KIN data frame
#' @param output.path Path to the directory where the output file will be saved.
#' @param output.id file name of the given network
#' @param functional_scores Computed functional scores if they were calculated (default: None)
#' @param beta FS filter threshold (default: 0.4)
#' @param gamma DIFF filter threshold (default: 1.0)
#' @return A data frame containing the correlation scores.
#' @examples
run_PCST <- function(
    baseline_KIN, DIFF_net, FS_net, DIFFandFS_net, output.path, 
    output.id, serine_threonine_kinase_data, tyrosine_kinase_data, gamma = 1.0) {
  
  # Helper function to compute PCST
  compute_PCST <- function(KIN, combination) {
    ## Preperation of dataframe for PCST
    # Adding helper edges from kinase substrates to kinases
    KIN$Type <- 'KS'
    helper_edges <- tibble(
      Source = unique(KIN$Target[which(KIN$Target_Uniprot %in% c(serine_threonine_kinase_data$kinase_name_mappings$`ACC#`, tyrosine_kinase_data$kinase_name_mappings$`ACC#`))]),
      Target = sapply(strsplit(Source, '_'), function(x) { x[1] }),
      Target_Uniprot = Target,
      f = gamma,
      Type = 'SK'
    )
    
    # Combine data and helper edges
    connected_KIN <- rbind(KIN, helper_edges)
    # reindex nodes and edges with integers
    node_names <- unique(c(connected_KIN$Source, connected_KIN$Target))
    nodes <- 0:(length(node_names)-1)
    names(nodes) <- node_names
    connected_KIN$Source_idx <- unname(sapply(connected_KIN$Source, function(x) { nodes[[x]]} ))
    connected_KIN$Target_idx <- unname(sapply(connected_KIN$Target, function(x) { nodes[[x]]} ))

    # create dir for PCST
    dir.create(paste0(output.path, '/PCST_networks/'), recursive = TRUE)

    # write node and edge file
    edge_file <- paste0(output.path, '/PCST_networks/', output.id, '_', combination, '_connectedNet.tsv')
    node_names_file <- paste0(output.path, '/PCST_networks/', output.id, '_', combination, '_connectedNet_nodeNames.tsv')
    write_tsv(
      connected_KIN, 
      edge_file
    )

    write_tsv(
      data.frame(
        name = names(nodes),
        index = unname(nodes),
        f = sapply(
          names(nodes), 
          function(x) {
            if (x %in% connected_KIN$Target) {
              unique(connected_KIN[which(connected_KIN$Target == x),]$f)
            } else {
              0
            }
          }
        )
      ),
      node_names_file
    )

    # Execute python script with PCST code
    system(
      paste0(
        'python -W ignore::DeprecationWarning src/python_src_scripts.py --net_file=', edge_file,
        ' --node_names_file=', node_names_file,
        ' --output_path=', output.path, '/PCST_networks/', output.id, '_', combination, '_net.tsv'
      )
    )
  }
  # without FS and DIFF filter
  compute_PCST(baseline_KIN[, c('Source', 'Target', 'Target_Uniprot', 'f')], 'PCST')
  compute_PCST(DIFF_net[, c('Source', 'Target', 'Target_Uniprot', 'f')], 'DIFFandPCST')
  if (!is.null(FS_net) & !is.null(DIFFandFS_net)) {
    compute_PCST(FS_net[, c('Source', 'Target', 'Target_Uniprot', 'f')], 'FSandPCST')
    compute_PCST(DIFFandFS_net[, c('Source', 'Target', 'Target_Uniprot', 'f')], 'DIFFandFSandPCST') 
  }
}

#' Run the pipeline to get the baseline KIN, functional scores, correlation filter and the PCST 
#' 
#' @param x1 Path to intensity data of condition 1 (default: NA).
#' @param x2 Path to intensity data of condition 2 (default: NA).
#' @param f Path to log2-transformed intensities (default: NA).
#' @param output_path Path to store the output (default: 'results').
#' @param output_id Output identifier (default: Key).
#' @param species Species name (default: Human).
#' @param paired_samples Bool variable that indicates if the samples are paired. Influences the computation of the log2-transformed intensities (default: FALSE).
#' @param log_intensities Bool variable that indicates if the intensities should be log tansformed if they are given as matrices (default: TRUE).
#' @param alpha Threshold for the baseline KIN percentile scoring (default: 0.9).
#' @param n Threshold for the number of edges in the baseline KIN (default: 15).
#' @param beta Threshold for the FS filter (default: 0.4).
#' @param gamma Threshold for the DIFF filter (default: 1.0).
#' @param delta Threshold for the CORR filter (default: 0.8).
#' @param epsilon Threshold for the significance test of the CORR filter (default: 0.05).
#' @param m Threshold for the number of samples needed for correlation calculation (default: 10).
#' @return A data frame containing the correlation scores.
#' @examples
#' # Example usage
#' run_pipeline('intensity_path', "./output/", disease_condition_NAME, control_condition_NAME)
run_KINference <- function(
  x1.path = NA, x0.path = NA, f.path = NA, output.path = 'results', output.id = 'key', species = 'Human', paired_samples = F,
  log_intensities = T, alpha = 0.9, n = 15, beta = 0.4, gamma = 1.0, delta = 0.8, epsilon = 0.05, m = 10 
  ) {
  
  message('Loading kinase data and preparing data!')
  # load kinase data
  serine_threonine_kinase_data <- load_kinase_data(
    kinase_motifs.path = './data/kinase_data/serine_threonine_kinases/kinase_motifs.csv',
    kinase_aprior_distributions.path = './data/kinase_data/serine_threonine_kinases/apriori_distributions/',
    kinase_name_mappings.path = './data/kinase_data/serine_threonine_kinases/kinase_name_mappings.tsv'
  )
  
  tyrosine_kinase_data <- load_kinase_data(
    kinase_motifs.path = './data/kinase_data/tyrosine_kinases/kinase_motifs.csv',
    kinase_aprior_distributions.path = './data/kinase_data/tyrosine_kinases/apriori_distributions/',
    kinase_name_mappings.path = './data/kinase_data/tyrosine_kinases/kinase_name_mappings.tsv'
  )
  
  # preprocessing data
  if (!is.na(x1.path) & !is.na(x0.path)) {
    
    x1 <- read_tsv(x1.path)
    x1 <- column_to_rownames(x1, 'Protein')
    x1 <- x1[, order(colnames(x1))]
    
    x0 <- read_tsv(x0.path)
    x0 <- column_to_rownames(x0, 'Protein')
    x0 <- x0[, order(colnames(x0))]
    
    if (paired_samples) {
      # adding NAs of one matrix to the other to enforce paired samples
      x1[which(is.na(x0), arr.ind = T)[1], which(is.na(x0), arr.ind = T)[2]]  <- NA
      x0[which(is.na(x1), arr.ind = T)[1], which(is.na(x1), arr.ind = T)[2]]  <- NA
    }
    # filter out all Proteins with no measurements (only NAs in the row)
    x1 <- x1[rowSums(is.na(x1)) != ncol(x1), ]
    x0 <- x0[rowSums(is.na(x0)) != ncol(x0), ]
    
    # filter out all proteins that are only measured in one condition
    x1 <- x1[which(rownames(x1) %in% rownames(x0)), ]
    x0 <- x0[which(rownames(x0) %in% rownames(x1)), ]
    
    # Order rownames just in case
    x1 <- x1[order(rownames(x1)),]
    x0 <- x0[order(rownames(x0)),]
    
    # log2 transform if needed
    if (!log_intensities) {
      x1 <- log2(x1)
      x0 <- log2(x0)
    }
    
    # Creating meta data
    md <- tibble(
      Protein = rownames(x1),
      Uniprot = sapply(strsplit(Protein, '_'), function(x) { x[1] }),
      AAPos = sapply(strsplit(Protein, '_'), function(x) { x[2] }),
      AA = sapply(strsplit(AAPos, '_'), function(x) { str_sub(x, 1, 1) }),
      Pos = sapply(strsplit(AAPos, '_'), function(x) { as.numeric(str_sub(x, 2)) }),
      Uniprot_Pos = paste0(Uniprot, '_', Pos)
    )
    
    if (species == 'Mus Musculus') {
      uniprot_MusMusculus <- md$Uniprot
      md$Uniprot <- orthology_translate_column(data = md, column = 'Uniprot', target_organism = 'human', source_organism = 'mouse')$Uniprot_9606
      md$Uniprot_musMusculus <- uniprot_MusMusculus
    }
    
    # remove all proteins that are not phosphorylated at serine (S) or threonine (T)
    ind_to_keep <- which(md$AA %in% c('S', 'T', 'Y'))
    x1 <- x1[ind_to_keep, ]
    x0 <- x0[ind_to_keep, ]
    md <- md[ind_to_keep, ]
    
    # compute log2FC
    if (!paired_samples) {
      md$f <- rowMeans(x1, na.rm = T) - rowMeans(x0, na.rm = T)
    } else {
      md$f <- rowMeans(x1 - x0, na.rm = T)
    }
  } else if(!is.na(f.path)) {
    f <- read_tsv(f.path)
    f <- column_to_rownames(f, 'Protein')
    
    md <- tibble(
      Protein = rownames(f),
      Uniprot = sapply(strsplit(Protein, '_'), function(x) { x[1] }),
      AAPos = sapply(strsplit(Protein, '_'), function(x) { x[2] }),
      AA = sapply(strsplit(AAPos, '_'), function(x) { str_sub(x, 1, 1) }),
      Pos = sapply(strsplit(AAPos, '_'), function(x) { as.numeric(str_sub(x, 2)) }),
      Uniprot_Pos = paste0(Uniprot, '_', Pos),
      f = f[[1]]
    )
    
    if (species == 'Mus Musculus') {
      uniprot_MusMusculus <- md$Uniprot
      md$Uniprot <- orthology_translate_column(data = md, column = 'Uniprot', target_organism = 'human', source_organism = 'mouse')$Uniprot_9606
      md$Uniprot_musMusculus <- uniprot_MusMusculus
    }
    
    md <- md[ which(md$AA %in% c('S', 'T', 'Y')), ]
    
  } else{
    stop('X1 and X0 or f have to be given!')
  }
  
  message('Inferring baseline KIN!')
  baseline_serine_threonine_KIN <- infer_baseline_serine_threonine_KIN(
    md, 
    output.path, 
    output.id, 
    serine_threonine_kinase_data
    )
  
  baseline_tyrosine_KIN <- infer_baseline_tyrosine_KIN(
    md, 
    output.path, 
    output.id, 
    tyrosine_kinase_data
  )
  
  message('Computing kinase enrichments!')
  # compute enriched kinases
  kinase_enrichment(
    md = md[which(md$AA %in% c('S', 'T')),],
    baseline_KIN = baseline_serine_threonine_KIN, 
    output.path = output.path, 
    output.id = paste0(output.id, '_serine_threonine_kinases'), 
    kinase_data = serine_threonine_kinase_data,
    gamma = gamma
    )
  
  kinase_enrichment(
    md = md[which(md$AA %in% c('Y')),],
    baseline_KIN = baseline_tyrosine_KIN, 
    output.path = output.path, 
    output.id = paste0(output.id, '_tyrosine_kinases'), 
    kinase_data = tyrosine_kinase_data,
    gamma = gamma
  )
  
  baseline_KIN <- rbind(baseline_serine_threonine_KIN, baseline_tyrosine_KIN)

  # Perform functional scoring (Only if species == HUMAN)
  if (species %in% c('Human', 'Mus Musculus')) {
    message('Computing FS!')
    functional_scores <- compute_FS_filter(
      md = md,
      output.path = output.path,
      output.id
      )
    
    functional_scores <- functional_scores[which(functional_scores$sites %in% md$Uniprot_Pos),]
    functional_scores <- functional_scores[order(functional_scores$sites), ]
    md <- md[order(md$Uniprot_Pos),]
    md$FS <- functional_scores$probabilities
    md <- md[order(md$Protein),]
    baseline_KIN$FS <- sapply(baseline_KIN$Target, function(x) { md$FS[md$Protein == x] })
  } else {
    functional_scores <- NA
  }

  dir.create(paste0(output.path, '/node_filtered_networks/'), recursive = TRUE)
  # Adding SK edges for the outputs
  add_SK_edges <- function(KIN) {
    KIN$Type = 'KS'
    KIN_SK <- tibble(
      Source =  unique(KIN$Target[which(KIN$Target_Uniprot %in% c(serine_threonine_kinase_data$kinase_name_mappings$`ACC#`, tyrosine_kinase_data$kinase_name_mappings$`ACC#`))]),
      Target = sapply(strsplit(Source, '_'), function(x) { x[1] }),
      Target_Uniprot = Target,
      ModifiedSequence = '',
      log2Score = 0,
      percentileScore = 0,
      percentileRank = 0,
      f = 0,
      FS = 0,
      Type = 'SK'
    )
    return(rbind(KIN, KIN_SK))
  }
  DIFF_net <- baseline_KIN[which(abs(baseline_KIN$f) >= gamma), ]
  write_tsv(
    add_SK_edges(DIFF_net),
    paste0(output.path, '/node_filtered_networks/', output.id, '_DIFFnet.tsv')
  )

  if (species %in% c('Human', 'Mus Musculus')) {
    FS_net <- baseline_KIN[which(baseline_KIN$Target %in% md$Protein[md$FS >= beta]), ]  
    DIFFandFS_net <- baseline_KIN[
      which(abs(baseline_KIN$f) >= gamma & baseline_KIN$Target %in% md$Protein[md$FS >= beta]),
    ]
    write_tsv(
      add_SK_edges(FS_net),
      paste0(output.path, '/node_filtered_networks/', output.id, '_FSnet.tsv')
    )
    write_tsv(
      add_SK_edges(DIFFandFS_net),
      paste0(output.path, '/node_filtered_networks/', output.id, '_DIFFandFSnet.tsv')
    )
  } else {
    FS_net <- NA
    DIFFandFS_net <- NA
  }
  
  # Perform CORR Filter
  if (!is.na(x1.path) & !is.na(x0.path)) {
    if (ncol(x1) >= m & ncol(x0) >= m) {
      message('Computing CORR!')
      # 1) For condition X1
      corr_x1 <- compute_CORR_filter(
        x = x1,
        md = md,
        baseline_KIN = baseline_KIN,
        output.path = output.path,
        output.id = paste0(output.id, '_X1'),
        serine_threonine_kinase_data,
        tyrosine_kinase_data,
        delta = delta,
        epsilon = epsilon,
        m = m  
      )
      corr_x1_DIFF_net <- corr_x1[which(paste0(corr_x1$Source, '_', corr_x1$Target) %in% paste0(DIFF_net$Source, '_', DIFF_net$Target)),]
      corr_x1_FS_net <- corr_x1[which(paste0(corr_x1$Source, '_', corr_x1$Target) %in% paste0(FS_net$Source, '_', FS_net$Target)),] 
      corr_x1_DIFFandFS_net <- corr_x1[which(paste0(corr_x1$Source, '_', corr_x1$Target) %in% paste0(DIFFandFS_net$Source, '_', DIFFandFS_net$Target)),] 
      write_tsv(
        add_SK_edges(corr_x1_DIFF_net),
        paste0(output.path, '/correlation_networks/', output.id, '_X1_DIFFnet.tsv')
      )
      write_tsv(
        add_SK_edges(corr_x1_FS_net),
        paste0(output.path, '/correlation_networks/', output.id, '_X1_FSnet.tsv')
      )
      write_tsv(
        add_SK_edges(corr_x1_DIFFandFS_net),
        paste0(output.path, '/correlation_networks/', output.id, '_X1_DIFFandFSnet.tsv')
      )
      
      # 2) For condition X0
      corr_x0 <- compute_CORR_filter(
        x = x0,
        md = md,
        baseline_KIN = baseline_KIN,
        output.path = output.path,
        output.id = paste0(output.id, '_X0'),
        serine_threonine_kinase_data,
        tyrosine_kinase_data,
        delta = delta,
        epsilon = epsilon,
        m = m    
      )
      corr_x0_DIFF_net <- corr_x0[which(paste0(corr_x0$Source, '_', corr_x0$Target) %in% paste0(DIFF_net$Source, '_', DIFF_net$Target)),]
      corr_x0_FS_net <- corr_x0[which(paste0(corr_x0$Source, '_', corr_x0$Target) %in% paste0(FS_net$Source, '_', FS_net$Target)),] 
      corr_x0_DIFFandFS_net <- corr_x0[which(paste0(corr_x0$Source, '_', corr_x0$Target) %in% paste0(DIFFandFS_net$Source, '_', DIFFandFS_net$Target)),]
      write_tsv(
        add_SK_edges(corr_x0_DIFF_net),
        paste0(output.path, '/correlation_networks/', output.id, '_X0_DIFFnet.tsv')
      )
      write_tsv(
        add_SK_edges(corr_x0_FS_net),
        paste0(output.path, '/correlation_networks/', output.id, '_X0_FSnet.tsv')
      )
      write_tsv(
        add_SK_edges(corr_x0_DIFFandFS_net),
        paste0(output.path, '/correlation_networks/', output.id, '_X0_DIFFandFSnet.tsv')
      )
    }
  }
  
  message('Computing PCST!')
  # Run PCST
  run_PCST(
    baseline_KIN = baseline_KIN,
    DIFF_net = DIFF_net,
    FS_net = FS_net,
    DIFFandFS_net = DIFFandFS_net,
    output.path = output.path,
    output.id = output.id,
    serine_threonine_kinase_data,
    tyrosine_kinase_data,
    gamma = gamma
  ) 
}