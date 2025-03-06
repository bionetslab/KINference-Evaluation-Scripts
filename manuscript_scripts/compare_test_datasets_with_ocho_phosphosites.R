library(tidyverse)

johnson_serine_threonine_substrates <- read_tsv('../data/kinase_data/serine_threonine_kinases/computing_aprior_distributions/substrates.tsv')
yaron_barir_tyrosine_substrates <- read_tsv('../data/kinase_data/tyrosine_kinases/computing_aprior_distributions/substrates.tsv')

bouhaddou_2023_x0 <- read_tsv('../data/Bouhaddou2023/intensities_VIC10h.tsv') 
bouhaddou_2023_x1 <- read_tsv('../data/Bouhaddou2023/intensities_Mock10h.tsv') 

print(paste0('Overlap between Johnson and Bouhaddou 2023 x0: ', bouhaddou_2023_x0 %>% filter((Protein %in% johnson_serine_threonine_substrates$Protein) | (Protein %in% yaron_barir_tyrosine_substrates$Protein)) %>% nrow()))
print(paste0('Overlap between Johnson and Bouhaddou 2023 x1: ', bouhaddou_2023_x1 %>% filter((Protein %in% johnson_serine_threonine_substrates$Protein) | (Protein %in% yaron_barir_tyrosine_substrates$Protein)) %>% nrow()))

wilkes_2015_f <- read_tsv('../data/Wilkes2015/intensities_log2fcTransformedVector.tsv')
print(paste0('Overlap between Johnson and Wilkes 2015 f: ', wilkes_2015_f %>% filter((Protein %in% johnson_serine_threonine_substrates$Protein) | (Protein %in% yaron_barir_tyrosine_substrates$Protein)) %>% nrow()))