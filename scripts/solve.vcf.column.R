library(readr)
library(purrr)


args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output <- args[2]
vcf_header <- read_lines(vcf_file) %>%
  purrr::keep(~grepl("^##", .))
vcf <- read_tsv(vcf_file, comment = "##") %>%
  dplyr::mutate(REF = "A", ALT = "T")
file.create(output)
write_lines(vcf_header, output, append = TRUE)
write_tsv(vcf, output, col_names = TRUE, append = TRUE)

