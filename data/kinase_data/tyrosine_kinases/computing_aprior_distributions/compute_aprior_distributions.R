library(tidyverse)
library(data.table)
library(UniProt.ws)

compute_apriori_distribution <- function(
    md, kinase_data) {

  kinases <- kinase_data$KINASE
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
    scores_log2 <- log2(apply(kinase_data[, ..aa], 1, prod))
    names(scores_log2) <- kinases
    scores_log2
  }
  scores_log2_list <- lapply(sequences_list, get_scores_log2, kinase_data = kinase_data)
  return(scores_log2_list)
}

substrates <- read_tsv('./data/kinase_data/tyrosine_kinases/computing_aprior_distributions/substrates.tsv')
kinase_data <- fread('./data/kinase_data/tyrosine_kinases/kinase_motifs.csv')
md <- tibble(
  Protein = substrates$Protein,
  Uniprot = sapply(strsplit(Protein, '_'), function(x) { x[1] }),
  AAPos = sapply(strsplit(Protein, '_'), function(x) { x[2] }),
  AA = sapply(strsplit(AAPos, '_'), function(x) { str_sub(x, 1, 1) }),
  Pos = sapply(strsplit(AAPos, '_'), function(x) { as.numeric(str_sub(x, 2)) }),
  Uniprot_Pos = paste0(Uniprot, '_', Pos)
)

log2_scores <- compute_apriori_distribution(md, kinase_data)

dir.create('./data/kinase_data/tyrosine_kinases/apriori_distributions/')

for (i in 1:nrow(kinase_data)){
  elements <- sapply(log2_scores, function(x) x[[i]])
  write_delim(data.frame(elements), paste0('./data/kinase_data/tyrosine_kinases/apriori_distributions/', kinase_data$KINASE[i], '.txt'), col_names = FALSE)
}

