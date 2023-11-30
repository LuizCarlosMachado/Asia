library(readr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output <- args[2]
vcf <- read_tsv(vcf_file, comment = "##") %>%
   dplyr::mutate(POS_Gen = POS/10^6) %>%
   dplyr::select(`#CHROM`, ID, POS_Gen, POS)
file.create(output)
write_tsv(vcf,output, col_names = F, append = TRUE)
