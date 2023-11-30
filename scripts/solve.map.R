library(readr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
map_file <- args[1]
output <- args[2]
map <- data.table::fread(map_file) %>%
   dplyr::mutate(V3 = V4/10^6)
file.create(output)
write_tsv(map, output, col_names = FALSE, append = TRUE)
