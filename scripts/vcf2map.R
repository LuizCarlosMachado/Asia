# Load necessary libraries
library(readr)  # Used for reading and writing tabular data
library(purrr)  # Functional programming tools, not explicitly used in this snippet

# Capture command line arguments for dynamic script usage
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]  # Input VCF file path
output <- args[2]    # Output file path

# Process the VCF file
vcf <- read_tsv(vcf_file, comment = "##") %>%
  dplyr::mutate(POS_Gen = POS/10^6) %>%  # Add a new column `POS_Gen` converting POS to megabases
  dplyr::select(`#CHROM`, ID, POS_Gen, POS)  # Select specific columns including the new `POS_Gen`

# Create the output file if it doesn't already exist
file.create(output)

# Write the processed data to the output file
# `col_names = F` indicates that column names will not be written to the file
# `append = TRUE` allows appending to the file if it already exists
write_tsv(vcf, output, col_names = FALSE, append = TRUE)
