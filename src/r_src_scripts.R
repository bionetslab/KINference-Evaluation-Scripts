library(data.table)
library(stringr)
library(dplyr)
library(foreach)
library(doParallel)
library(UniProt.ws)
library(funscoR)
library(knitr)
library(OmnipathR)


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
load_kinase_data <- function(kinase_motifs.path = './data/kinase_data/kinase_motifs.csv',
                             kinase_aprior_distributions.path = './data/kinase_data/apriori_distributions/',
                             kinase_scaling_factors.path = './data/kinase_data/substrate_scoring_scaling_factors.csv',
                             kinase_name_mappings.path = './data/kinase_data/kinase_name_mappings.tsv',
                             kinase_name_mappings_mouse.path = './data/kinase_data/kinase_name_mappings_mouse.tsv') {
  
  # reading in kinase motif data
  kinase_motifs <- fread(kinase_motifs.path)

  # reading in mapping of kinases to uniprot id file
  kinase_names_to_uniprotmapping <- fread(kinase_name_mappings.path)
  kinase_names_to_uniprotmapping <- kinase_names_to_uniprotmapping[which(kinase_names_to_uniprotmapping$PROTEIN %in% kinase_motifs$KINASE),]
  
  kinase_names_to_uniprotmapping_mouse <- fread(kinase_name_mappings_mouse.path)
  kinases <- kinase_names_to_uniprotmapping[['ACC#']]

  # reading in kinase scaling factors (as used in the paper)
  if (!is.na(kinase_scaling_factors.path)) {
    kinase_scaling_factors <- fread(kinase_scaling_factors.path)$scale
  }

  # Computing a-priori distributions of kinases based on log2 score database
  kinase_aprior_distributions <- list()
  for (i in seq_len(nrow(kinase_motifs))) {
    row <- kinase_motifs[i, ]
    kin <- row$KINASE
    kin_aprior_file_name <- paste(kinase_aprior_distributions.path, kin, '.txt', sep = '')
    aprior_distribution <- ecdf(fread(kin_aprior_file_name)[[1]])
    kinase_aprior_distributions <- append(kinase_aprior_distributions, aprior_distribution)
  }
  names(kinase_aprior_distributions) <- kinases
  
  return(list(
    kinase_motifs = kinase_motifs, kinase_aprior_distributions = kinase_aprior_distributions, 
    kinase_scaling_factors = kinase_scaling_factors, kinase_name_mappings = kinase_names_to_uniprotmapping, kinase_name_mappings_mouse = kinase_names_to_uniprotmapping_mouse
  ))

}

#' Perform substrate scoring
#' 
#' This function calculates substrate scores as described in https://doi.org/10.1038/s41586-022-05575-3 Fig. 3
#' 
#' @param intensities A data frame containing intensities.
#' @param output.path Path to the directory where the output file will be saved.
#' @param output_filename Name of the output file.
#' @param kinase_data Kinase information loaded by the function load_kinase_data()
#' @param percentile_rank_threshold Threshold for percentile rank (default: 15).
#' @param alpha Threshold for percentile scores based of the percentage (default: 0.9).
#' @return A data frame containing the calculated substrate scores.
#' @examples
#' # Example usage
#' substrate_scoring(intensities_df, "./output/", "cond1_vs_cond2.tsv", kinase_motifs_df, kinase_apriori_distributions_list, kinase_scaling_factors_df, 15)
substrate_scoring <- function(
    intensities, output.path, output_filename, kinase_data, 
    percentile_rank_threshold = 15, alpha = 0.9) {
  
  kinases <- kinase_data$kinase_name_mappings$`ACC#`
  
  # Loading sequences from Uniprot  
  sequences <- uniprotSequencesFromWeb(unique(intensities$LeadingProtein))
  
  # Mapping sequences to conditions
  intensities$Seq <- sapply(intensities$LeadingProtein, function(x) { unlist(sequences[which(sequences$Entry == x), 'Sequence']) })
  intensities <- intensities[which(nchar(intensities$Seq) > 0),]
  # Removing all entries where the given phosphorylation position is not a 'S' or 'T'
  # Kinase motif data is only for Threonine and Seronine phosphorylation sites
  intensities <- intensities[which(str_sub(intensities$Seq, intensities$Pos, intensities$Pos) %in% c('S', 'T')),]
  intensities <- intensities[which(!grepl(';', intensities$Protein)),]
  intensities <- intensities[which(!is.infinite(intensities$log2FC)),]
  # create output directory for substrate scores  
  dir.create(paste0(output.path, '/substrate_scores/'), recursive = TRUE)

  # start parallel cluster for faster computation of substrate scores
  cores <- detectCores()
  cores <- min(9, cores[1])
  cl <- makeCluster(cores - 1) # not to overload your computer
  registerDoParallel(cl)

  # Substrate scoring
  sub_scores = foreach(i=1:nrow(intensities), .combine=rbind, .packages = c('data.table', 'stringr', 'dplyr')) %dopar% {
    # for (i in 1:nrow(intensities)) { # can be used for debugging
    row_i <- intensities[i, ]
    target <- row_i$Protein
    seq <- row_i$Seq
    pos <- row_i$Pos
    if (pos > nchar(seq)) {
      return(NULL)
    }
    # Get sequence around the phosphorylation site (-5 and +4 positions)
    if ((pos > 5) && ((pos + 4) < nchar(seq))) {
      # not at the borders of the sequence
      seq <- str_sub(seq, pos - 5, pos + 4)
      zero_idx <- 6
      seq_split <- str_extract_all(seq, stringr::boundary('character'))[[1]]
      indices <- seq(-5, 4, 1)
    } else if (pos > 5) {
      # at the right border of the sequence
      seq <- str_sub(seq, pos - 5)
      zero_idx <- 6
      seq_split <- str_extract_all(seq, stringr::boundary('character'))[[1]]
      indices <- seq(-5, nchar(seq)-6, 1)
    } else if ((pos + 4) < nchar(seq)) {
      # at the left border of the sequence
      seq <- str_sub(seq, 1, pos + 4)
      zero_idx <- (nchar(seq) - 4)
      seq_split <- str_extract_all(seq, stringr::boundary('character'))[[1]]
      indices <- seq(-(nchar(seq) - 5), 4, 1)
    } else {
      # weird short sequence -> skip it
      return(NULL) # use 'next' when debugging
    }

    # log2 of scaled product of matching columns in PSSM matrix to sequence
    aa <- paste0(indices, seq_split)[-zero_idx]
    scores_log2 <- log2(apply(kinase_data$kinase_motifs[, ..aa], 1, prod) / kinase_data$kinase_scaling_factors)
    
    # compute percentile scores
    names(scores_log2) <- kinases
    percentile_scores <- sapply(names(kinase_data$kinase_aprior_distributions), function(x) { kinase_data$kinase_aprior_distributions[[x]](scores_log2[[x]]) })
    
    # append new rows to dataframe
    percentile_scores <- percentile_scores[order(-percentile_scores)]
    
    edges_threshold <- min(percentile_rank_threshold, length(which(percentile_scores >= alpha)))
    # Check if any edge satisfies the thresholds
    if(edges_threshold == 0) {
      return(NULL)
    }
    
    new_rows <- data.frame(
      Source = names(percentile_scores)[1:edges_threshold], 
      Target = rep(target, edges_threshold), 
      Target_LeadingProtein = rep(row_i$LeadingProtein, edges_threshold),
      ModifiedSequence = rep(seq, edges_threshold), 
      log2Score = scores_log2[names(percentile_scores)[1:edges_threshold]], 
      percentileScore = percentile_scores[1:edges_threshold], 
      percentileRank = seq(1,edges_threshold,1), 
      log2FC = rep(row_i$log2FC, edges_threshold)
    )
    new_rows
  }

  # write the substrate scoring file to the results path
  fwrite(sub_scores, paste0(output.path, '/substrate_scores/', output_filename, '.tsv'), sep = '\t')

  # stop the parallel cluster
  stopCluster(cl)

  return(sub_scores)

}


#' Perform kinase enrichment analysis
#' 
#' This function performs kinase enrichment analysis based on intensities and substrate scores,
#' and writes the results to a specified output file.
#' 
#' @param intensities A data frame containing intensities.
#' @param sub_scores A data frame containing substrate scores.
#' @param output.path Path to the directory where the output file will be saved.
#' @param output_filename Name of the output file.
#' @param gamma Threshold for log2 fold change.
#' @param kinase_data Kinase information loaded by the function load_kinase_data()
#' @return A data frame containing the results of the kinase enrichment analysis.
#' @examples
#' # Example usage
#' kinase_enrichment(intensities_df, substrate_scores_df, "./output/", "cond1_vs_cond2.tsv", 1.0)
kinase_enrichment <- function(
    intensities, sub_scores, output.path, 
    output_filename, kinase_data, gamma = 1.0) {
  
  kinases <- kinase_data$kinase_name_mappings$`ACC#`
  # creating output directory for enrichments
  dir.create(paste0(output.path, '/enrichments/'), recursive = TRUE)
  # Upregulated set := Phosphosites with log2FC > gamma
  # Downregulated set := Phosphosites with log2FC < - (gamma)
  # Background set := Phosphosites with -(gamma) <= log2FC <= gamma
  upregulated_intensities <- intensities[which(intensities$log2FC > gamma),]
  downregulated_intensities <- intensities[which(intensities$log2FC < gamma),]
  background_intensities <- intensities[which(abs(intensities$log2FC) <= gamma),]

  # Getting all interactions in each set 
  upregulated_interactions <- sub_scores[which(sub_scores$Target %in% upregulated_intensities$Protein),]
  downregulated_interactions <- sub_scores[which(sub_scores$Target %in% downregulated_intensities$Protein),]
  background_interactions <- sub_scores[which(sub_scores$Target %in% background_intensities$Protein),]

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

  up_set_size <- nrow(upregulated_intensities)
  down_set_size <- nrow(downregulated_intensities)
  background_set_size <- nrow(background_intensities)
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
                      protein = sapply(up_counts$Var1, function(x) {kinase_data$kinase_name_mappings$PROTEIN[which(kinase_data$kinase_name_mappings$`ACC#` == x)]}),
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

  fwrite(result, paste0(output.path, '/enrichments/', output_filename, '.tsv'), sep = '\t')
  return(result)
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
functional_scoring <- function(intensities, output.path) {

  data_targets <- data.frame(
    acc=intensities$LeadingProtein,
    position=intensities$Pos,
    residue=intensities$AA
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
  
  fwrite(all_scores, paste0(output.path, '/functionalScores.tsv'), sep='\t')
  
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
#' @param output_filename Name of the output file.
#' @param m Threshold for the number of samples needed for correlation calculation (default: 10).
#' @return A data frame containing the correlation scores.
#' @examples
#' # Example usage
#' correlation_scoring(intensities_df, substrate_scores_df, "./output/", "correlation_scores.csv", 10)
correlation_scoring <- function(
  intensities, sub_scores, condition, 
  output.path, m = 10) {
  
  dir.create(paste0(output.path, '/correlation_networks/'), recursive = TRUE)
  intensities <- intensities[which(intensities$Condition == condition),]

  start <- Sys.time()
  # only keeping the edges that have an intensity measure on the source (the inferred network from the substrate scoring is kinase -> PhosphoSite not PhosphoSite -> PhosphoSite)
  sub_scores <- sub_scores[sub_scores$Source %in% intensities$LeadingProtein, ]
  
  cores <- detectCores()
  cores <- min(9, cores[1])
  cl <- makeCluster(cores-1) #not to overload your computer
  registerDoParallel(cl)
  
  # Computing correlation score and differential correlation score
  correlation_network = foreach(i=1:nrow(sub_scores), .combine = rbind, .packages = c('data.table', 'stringr', 'dplyr')) %dopar% {
    df <- data.frame()
    row <- sub_scores[i, ]
    target <- row$Target
    source <- row$Source
    
    target_intensities <- intensities[which(intensities$Protein == sub_scores[i, ]$Target), ]
    target_subjects <- target_intensities$Subject
    
    sources <- unique(intensities[which(intensities$LeadingProtein == source),]$Protein)
    for (j in seq_len(length(sources))) {
      
      source_intensities <- intensities[which(intensities$Protein == sources[j]), ] 
      source_subjects <- source_intensities$Subject
      
      subjects_measured_in_both <- intersect(source_subjects, target_subjects)
      
      if (length(subjects_measured_in_both) >= m) {
        
        x <- source_intensities[which(source_subjects %in% subjects_measured_in_both),]
        y <- target_intensities[which(target_subjects %in% subjects_measured_in_both),]
        
        x <- x[order(x$Subject),]$LogIntensity
        y <- y[order(y$Subject),]$LogIntensity
        
        data <- data.frame(x = x, y = y)
        correlation_test <- cor.test(data$x, data$y, method = "pearson")
        if (!is.na(correlation_test[["p.value"]]) && correlation_test[["p.value"]] <= 0.05) {
          corCoef <- unname(correlation_test[['estimate']])
          new_rows <- data.frame(
            Source = sources[j],
            Target = target,
            Target_LeadingProtein = row$Target_LeadingProtein,
            log2Score = row$log2Score,
            percentileScore = row$percentileScore,
            percentileRank = row$percentileRank,
            log2FC = row$log2FC,
            corCoef = corCoef
          )
          
          df <- rbind(df, new_rows)
        }
      } 
    }
    df
  }
  
  stopCluster(cl)
  
  end <- Sys.time()
  

  fwrite(correlation_network, paste0(output.path, '/correlation_networks/', condition, '.tsv'), sep = "\t")
  message(paste0('Finished correlation scoring for ', condition, ' Time: ', end - start))
  
}

#' Get intensities for pairwise comparison of conditions
#' 
#' This function takes a data frame of intensities and two condition labels (c1 and c2) and
#' returns the intensities corresponding to these two conditions for pairwise comparison with log2FC = c1$logIntensity - c2$logIntensity.
#' 
#' @param intensities A data frame containing intensities.
#' @param disease_condition Label for the disease condition.
#' @param control_condition Label for the control condition.
#' @return A list containing the intensities for the two conditions.
#' @examples
#' # Example usage
#' get_intensities_for_pairwise_comparison(intensities_df, "condition1", "condition2")
get_intensities_for_pairwise_comparison <- function(intensities, disease_condition, control_condition) {

  # only keeping entries for conditions c1 and c2
  intensities_conditions <- intensities[which(intensities$Condition %in% c(disease_condition, control_condition)),]
  
  # Computing mean log intensities for condition c1 and c2
  intensities_diseaseCondition <- intensities_conditions[intensities_conditions$Condition == disease_condition,] %>%
    dplyr::group_by(Protein, ProteinPhoSite, LeadingProtein, AA, Pos) %>%
    dplyr::summarise(MeanLogIntensity = mean(LogIntensity, na.rm = 'TRUE')) %>%
    as.data.frame()

  intensities_controlCondition <- intensities_conditions[intensities_conditions$Condition == control_condition,] %>%
    dplyr::group_by(Protein, ProteinPhoSite, LeadingProtein, AA, Pos) %>%
    dplyr::summarise(MeanLogIntensity = mean(LogIntensity, na.rm = 'TRUE'), .groups = 'keep') %>%
    as.data.frame()

  # Can only compute log2FC between phospho sites that are measured in both conditions
  proteins_intersection <- intersect(intensities_diseaseCondition$Protein, intensities_controlCondition$Protein)
  intensities_diseaseCondition <- intensities_diseaseCondition[which(intensities_diseaseCondition$Protein %in% proteins_intersection),]
  intensities_controlCondition <- intensities_controlCondition[which(intensities_controlCondition$Protein %in% proteins_intersection),]

  # Sort them by ProteinPhoSite (remove any missmatches)
  intensities_diseaseCondition <- intensities_diseaseCondition[order(intensities_diseaseCondition$Protein), ]
  intensities_controlCondition <- intensities_controlCondition[order(intensities_controlCondition$Protein), ]

  # Compute mean log2FC
  intensities_diseaseCondition$log2FC <- intensities_diseaseCondition$MeanLogIntensity - intensities_controlCondition$MeanLogIntensity

  intensities_diseaseCondition <- intensities_diseaseCondition[!is.infinite(intensities_diseaseCondition$log2FC),]

  return(intensities_diseaseCondition)
}

#' Run the PCST algorithm on the given baseline KIN
#' 
#' @param sub_scores Baseline KIN data frame
#' @param output.path Path to the directory where the output file will be saved.
#' @param output_filename file name of the given network
#' @param functional_scores Computed functional scores if they were calculated (default: None)
#' @param beta FS filter threshold (default: 0.4)
#' @param gamma DIFF filter threshold (default: 1.0)
#' @return A data frame containing the correlation scores.
#' @examples
#' # Example usage
#' run_pipeline(sub_scores, "./output/", output_filename)
run_PCST <- function(sub_scores, output.path, output_filename, kinase_data, functional_scores, beta = 0.4, gamma = 1.0) {
  
  # Helper function to compute PCST
  compute_PCST <- function(sub_score_net, file_name) {
    ## Preperation of dataframe for PCST
    sub_score_net$edge_type <- 'data'
    # Adding helper edges from kinase substrates to kinases
    connected_sub_score_net <- data.frame()
    for (i in 1:nrow(sub_score_net)) {
      if(sub_score_net$Target_LeadingProtein[i] %in% kinase_data$kinase_name_mappings$`ACC#`) {
        connected_sub_score_net <- rbind(
          connected_sub_score_net, 
          data.frame(
            Source = sub_score_net$Target[i], 
            Target = sub_score_net$Target_LeadingProtein[i],
            Target_LeadingProtein = sub_score_net$Target_LeadingProtein[i],
            log2FC = 0,
            edge_type = 'helper'
          )
        )
      }
    }

    # Combine data and helper edges
    connected_sub_score_net <- rbind(sub_score_net, unique(connected_sub_score_net))
    # reindex nodes and edges with integers
    node_names <- unique(c(connected_sub_score_net$Source, connected_sub_score_net$Target))
    nodes <- 0:(length(node_names)-1)
    names(nodes) <- node_names
    connected_sub_score_net$Source_idx <- unname(sapply(connected_sub_score_net$Source, function(x) { nodes[[x]]} ))
    connected_sub_score_net$Target_idx <- unname(sapply(connected_sub_score_net$Target, function(x) { nodes[[x]]} ))

    # create dir for PCST
    dir.create(paste0(output.path, '/PCST_networks/'), recursive = TRUE)

    # write node and edge file
    edge_file <- paste0(output.path, '/PCST_networks/', file_name, '_connected_subscore_net.tsv')
    node_names_file <- paste0(output.path, '/PCST_networks/', file_name, '_connected_subscore_net_nodeNames.tsv')
    fwrite(
      connected_sub_score_net, 
      paste0(output.path, '/PCST_networks/', file_name, '_connected_subscore_net.tsv'),
      sep = '\t'
    )

    fwrite(
      data.frame(
        name = names(nodes),
        index = unname(nodes),
        log2FC = sapply(
          names(nodes), 
          function(x) {
            if (x %in% connected_sub_score_net$Target) {
              unique(connected_sub_score_net[which(connected_sub_score_net$Target == x),]$log2FC)
            } else {
              0
            }
          }
        )
      ),
      paste0(output.path, '/PCST_networks/', file_name, '_connected_subscore_net_nodeNames.tsv'),
      sep = '\t'
    )

    # Execute python script with PCST code
    system(
      paste0(
        'python -W ignore::DeprecationWarning src/python_src_scripts.py --net_file=', edge_file,
        ' --node_names_file=', node_names_file,
        ' --output_path=', output.path, '/PCST_networks/', file_name, '_PCST_network.tsv'
      )
    )
  }
  # without FS and DIFF filter
  sub_score_net <- sub_scores[, c('Source', 'Target', 'Target_LeadingProtein', 'log2FC')]
  compute_PCST(sub_score_net, output_filename)
  
  # with DIFF filter
  sub_score_net_withDIFFFilter <- sub_score_net[which(abs(sub_score_net$log2FC) >= gamma), ]
  compute_PCST(
    sub_score_net_withDIFFFilter, 
    paste0(output_filename, '_withDIFFfilter')
  )
  
  if (any(!is.na(functional_scores))) {
    # with FS filter
    sub_score_net_withFSFilter <- sub_score_net
    sub_score_net_withFSFilter$TargetPos <- sapply(strsplit(sub_score_net_withFSFilter$Target, '_'), function(x) { str_sub(x[[2]], 2) })
    sub_score_net_withFSFilter$Target_w_Pos <- paste0(sub_score_net_withFSFilter$Target_LeadingProtein, '_', sub_score_net_withFSFilter$TargetPos)
    probabilites <- functional_scores$probabilities
    sub_score_net_withFSFilter$functionalScore <- sapply(
      sub_score_net_withFSFilter$Target_w_Pos, function(x) {
        probabilites[which(functional_scores$sites == x)]
      }
    )
    
    compute_PCST(
      sub_score_net[which(sub_score_net_withFSFilter$functionalScore >= beta), ], 
      paste0(output_filename, '_withFSfilter')
    )
    
    # with DIFF and FS filter
    compute_PCST(
      sub_score_net[which(sub_score_net_withFSFilter$functionalScore >= beta & sub_score_net_withFSFilter$log2FC >= gamma), ], 
      paste0(output_filename, '_withFSandDIFFfilter')
    )
  }

}


#' Run the pipeline to get the baseline KIN, functional scores, correlation filter and the PCST 
#' 
#' @param intensity_path Path to intensitiy data.
#' @param output.path Path to the directory where the output file will be saved.
#' @param disease_condition Disease condition name
#' @param control_condition Control condition name
#' @param species Species name (default: HUMAN)
#' @param compute_CORR_filter Bool variable that indicates if the CORR filter should be calculated
#' @param m Threshold for the number of samples needed for correlation calculation (default: 10).
#' @param alpha Threshold for the baseline KIN percentile scoring (default: 0.9).
#' @param beta Threshold for the FS filter (default: 0.4).
#' @param gamma Threshold for the DIFF filter (default: 1.0) 
#' @return A data frame containing the correlation scores.
#' @examples
#' # Example usage
#' run_pipeline('intensity_path', "./output/", disease_condition_NAME, control_condition_NAME)
run_pipeline <- function(
  intensity.path, output.path, disease_condition, control_condition, 
  species = 'HUMAN', compute_CORR_filter = FALSE, m = 10, alpha = 0.9, beta = 0.4, gamma = 1.0) {
  
  message('Loading kinase data and preparing data.')
  # load kinase data
  kinase_data <- load_kinase_data()

  # define output_filename
  output_filename <- paste0(disease_condition, '_vs_', control_condition)

  # load intensities
  intensities_unproc <- fread(intensity.path)
  # stop if the necessary columns are not present
  if (!all(c('Protein', 'LogIntensity', 'Condition', 'Subject') %in% colnames(intensities_unproc))) {
    stop('Some of the columns "Protein", "LogIntensity", "Condition" or "Subject" were not found in the input file.\n These columns have to exist for downstream analysis.\n Refer to the example_run.R script or the README in the github for more information!')
  }

  # add convenience columns to the intensity df
  intensities_unproc$ProteinPhoSite <- sapply(strsplit(intensities_unproc$Protein, ';'), function(x) { x[1] })
  intensities_unproc$AA <- sapply(strsplit(intensities_unproc$ProteinPhoSite, '_'), function(x) { str_sub(x[2], 1, 1) })
  intensities_unproc$Pos <- sapply(strsplit(intensities_unproc$ProteinPhoSite, '_'), function(x) { as.numeric(str_sub(x[2], 2)) })
  intensities_unproc$LeadingProtein <- sapply(strsplit(intensities_unproc$ProteinPhoSite, '_'), function(x) { x[1] })

  # get intensities for the two conditions
  intensities <- get_intensities_for_pairwise_comparison(intensities_unproc, disease_condition, control_condition)

  message('Inferring baseline KIN.')
  # perform substrate scoring <-> Compute baseline KIN
  sub_scores <- substrate_scoring(
    intensities = intensities, 
    output.path = output.path, 
    output_filename = output_filename, 
    kinase_data =  kinase_data,
    alpha = alpha
  )

  message('Computing kinase enrichments (Johnson et al. 2023).')
  # compute enriched kinases
  kinase_enrichment(
    intensities = intensities, 
    sub_scores = sub_scores, 
    output.path = output.path, 
    output_filename = output_filename, 
    kinase_data = kinase_data,
    gamma = gamma
  )

  # Perform functional scoring (Only if species == HUMAN)
  if (species == 'HUMAN') {
    message('Computing functional scores.')
    functional_scores <- functional_scoring(
      intensities = intensities,
      output.path = output.path
    )
  } else {
    functional_scores <- NA
  }

  # Perform CORR Filter
  if (compute_CORR_filter) {
    message('Computing CORR.')
    # 1) For the disease condition
    correlation_scoring(
      intensities = intensities_unproc,
      sub_scores = sub_scores,
      condition = disease_condition,
      output.path = output.path,
      m = m    
    )
    # 2) For the control condition
    correlation_scoring(
      intensities = intensities_unproc,
      sub_scores = sub_scores,
      condition = control_condition,
      output.path = output.path,
      m = m    
    )
  }

  message('Computing PCST.')
  # Run PCST
  run_PCST(
    sub_scores = sub_scores, 
    output.path = output.path,
    output_filename = output_filename,
    kinase_data = kinase_data,
    functional_scores = functional_scores,
    beta = beta,
    gamma = gamma
  ) 
  
}

run_pipeline_withGivenLog2FC <- function(
    intensity.path, output.path, disease_condition, control_condition = 'Mock', 
    species = 'HUMAN', m = 10, alpha = 0.9, beta = 0.4, gamma = 1.0) {

  message('Loading kinase data and preparing data.')
  # load kinase data
  kinase_data <- load_kinase_data()
  
  # define output_filename
  output_filename <- paste0(disease_condition, '_vs_', control_condition)
  
  # load intensities
  intensities <- fread(intensity.path)
  # stop if the necessary columns are not present
  if (!all(c('Protein', disease_condition) %in% colnames(intensities))) {
    stop(paste0('One of the columns "Protein" or ', disease_condition, ' were not found in the input file.\n These columns have to exist for downstream analysis.\n Refer to the example_run.R script or the README in the github for more information!'))
  }  
  intensities <- intensities[, .SD, .SDcols = c('Protein', disease_condition)]
  colnames(intensities) <- c('Protein', 'log2FC')
  
  # add convenience columns to the intensity df
  intensities$ProteinPhoSite <- sapply(strsplit(intensities$Protein, ';'), function(x) { x[1] })
  intensities$AA <- sapply(strsplit(intensities$ProteinPhoSite, '_'), function(x) { str_sub(x[2], 1, 1) })
  intensities$Pos <- sapply(strsplit(intensities$ProteinPhoSite, '_'), function(x) { as.numeric(str_sub(x[2], 2)) })
  intensities$LeadingProtein <- sapply(strsplit(intensities$ProteinPhoSite, '_'), function(x) { x[1] })
  
  message('Inferring baseline KIN.')
  # perform substrate scoring <-> Compute baseline KIN
  sub_scores <- substrate_scoring(
    intensities = intensities, 
    output.path = output.path, 
    output_filename = output_filename, 
    kinase_data =  kinase_data,
    alpha = alpha
  )
  
  message('Computing kinase enrichments (Johnson et al. 2023).')
  # compute enriched kinases
  kinase_enrichment(
    intensities = intensities, 
    sub_scores = sub_scores, 
    output.path = output.path, 
    output_filename = output_filename, 
    kinase_data = kinase_data,
    gamma = gamma
  )
  

  # Perform functional scoring (Only if species == HUMAN)
  if (species == 'HUMAN') {
    message('Computing functional scores.')
    functional_scores <- functional_scoring(
      intensities = intensities,
      output.path = output.path
    )
  } else {
    functional_scores <- NA
  }
  
  # Skipping CORR as there are no sample specific log intensity measurements
  
  message('Computing PCST.')
  # Run PCST
  run_PCST(
    sub_scores = sub_scores, 
    output.path = output.path,
    output_filename = output_filename,
    kinase_data = kinase_data,
    functional_scores = functional_scores,
    beta = beta,
    gamma = gamma
  )

}
