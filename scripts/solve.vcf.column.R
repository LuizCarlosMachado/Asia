# Load the necessary libraries
library(readr)  # For reading and writing files
library(purrr)  # For functional programming tools

# Parse command line arguments to get the input VCF file and output file paths
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]  # The path to the input VCF file
output <- args[2]    # The path to the output file

# Extract header lines from the VCF file
vcf_header <- read_lines(vcf_file) %>%
  purrr::keep(~grepl("^##", .))  # Keep lines that start with "##", which are header comments in VCF files

# Read the VCF file, skipping header comments, and modify the REF and ALT columns
vcf <- read_tsv(vcf_file, comment = "##") %>%
  dplyr::mutate(REF = "A", ALT = "T")  # Replace all REF values with "A" and all ALT values with "T"

# Create the output file if it does not exist
file.create(output)

# Write the extracted header comments to the output file
write_lines(vcf_header, output, append = TRUE)

# Write the modified VCF data to the output file, including column names and appending to the file
write_tsv(vcf, output, col_names = TRUE, append = TRUE)

