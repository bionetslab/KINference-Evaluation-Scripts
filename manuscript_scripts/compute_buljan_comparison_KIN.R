library(tidyverse)
library(readxl)
library(Kinference)
# Load the data
data <- read_excel('../data/competitors/buljan2020_interactions.xlsx', na = c('NA', ''))
colnames(data) <- data[2, ]
data <- data %>%
  filter(!row_number() %in% c(1, 2)) %>%
  filter(!is.na(wdscore), !is.na(gfpratio)) %>%
  mutate(wdscore = as.numeric(wdscore), gfpratio = as.numeric(gfpratio)) %>%
  filter(wdscore >= 73.6, gfpratio >= 18.4)

# Get targeted proteins
targets_uniprotIDs <- data %>% select(Protein_id) %>% distinct() %>% pull()

# Get sequences
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

uniprotSequences <- uniprotSequencesFromWeb(targets_uniprotIDs)

# Get all possible phosphorylation positions
Kinference_input <- data.frame()
AAs <- c('S', 'T', 'Y')
for (i in 1:length(uniprotSequences$query)) {
    for (AA in AAs) {
        matches <- gregexpr(AA, uniprotSequences$Sequence[i])
        positions <- as.integer(matches[[1]])
        if (positions[1] != -1) {
            for (pos in positions) {
                Kinference_input <- rbind(
                    Kinference_input,
                    list(
                        Protein=paste0(uniprotSequences$query[i], '_', AA, pos),
                        f=1
                    )
                )
            }
        }
    }
}

write_tsv(Kinference_input, './buljan2020_baselineKIN_comparison_data.tsv')

# Run Kinference
Kinference::run_KINference(
    f.path='./buljan2020_baselineKIN_comparison_data.tsv',
    gamma=1.0,
    output.path='../results/baselineKIN_comparisons',
    output.id='buljan2020_comparison',
    apply.PCST=FALSE,
    apply.FS=TRUE,
    apply.CORR=FALSE,
    apply.DIFF=FALSE,
    species='Homo sapiens',
    custom_serine_threonine_kinase_data.path=NULL,
    custom_tyrosine_kinase_data.path=NULL,
    n=300
)
