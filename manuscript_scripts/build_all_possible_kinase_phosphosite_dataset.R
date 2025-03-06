library(dplyr)

#' Query UniProt Web Service for Protein Sequences
#'
#' This function queries the UniProt web service for protein sequences based on 
#' a vector of UniProt accession numbers.
#'
#' @param uniprots A character vector of UniProt accession numbers to query.
#' @param fields A character vector of fields to retrieve from the UniProt web service. 
#'        Default is c("accession", "sequence").
#'
#' @return A data.frame containing the queried information from UniProt, with an 
#'         additional column 'query' that matches the input accession numbers.
#' @noRd
.query_uniprot_ws_for_sequence <- function(
    uniprots, 
    fields = c("accession","sequence")
) {
    queries <- paste0("accession:", uniprots)

    res <- UniProt.ws::queryUniProt(queries, fields = fields)
    res <- res %>%
        dplyr::mutate(query = Entry) %>%
        dplyr::right_join(data.frame(Entry = uniprots), by = "Entry")
    return(res)
}

#' Retrieve UniProt Sequences from Web
#'
#' This function retrieves protein sequences from the UniProt database for a given set of UniProt IDs.
#'
#' @param uniprotIDs A character vector of UniProt IDs for which sequences are to be retrieved.
#' @param chunkSize An integer specifying the number of UniProt IDs to query in each batch. Default is 25.
#' @param fields A character vector specifying the fields to retrieve from UniProt. Default is c("accession", "sequence").
#'
#' @return A data.frame containing the retrieved information for the specified UniProt IDs.
#' @noRd
uniprotSequencesFromWeb <- function(
    uniprotIDs, 
    chunkSize = 25, 
    fields = c("accession","sequence")
    )
{
    uniqueUniprots <- unique(uniprotIDs)
    chunks <- split(uniqueUniprots, (1:length(uniqueUniprots))%/%chunkSize)
    infoMapList <- pbapply::pblapply(chunks, .query_uniprot_ws_for_sequence, fields = fields)
    infoMap <- dplyr::bind_rows(infoMapList)
    return (infoMap)
}

# Building kinase dummy database
kinases <- c(Kinference::serine_threonine_kinase_data$kinase_name_mappings[['ACC#']], Kinference::tyrosine_kinase_data$kinase_name_mappings[['ACC#']])
kinase_sequences <- uniprotSequencesFromWeb(kinases)
data <- data.frame()

for (kinase in kinases) {
    for (aa in c('S', 'T', 'Y')) {
        positions <- gregexpr(aa, kinase_sequences$Sequence[which(kinase_sequences$query == kinase)])[[1]]
        if (positions[1] != -1) {
            for (pos in positions) {
                data <- rbind(data, data.frame(
                    kinase = paste0(kinase, '_', aa, pos),
                    f = 1.0
                ))
            }
        }
    }
}

readr::write_tsv(data, './kinase_phosphosite_dataset.tsv')