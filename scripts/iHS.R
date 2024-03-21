# Load the necessary libraries
library(rehh)
library(dplyr)
library(readr)

# Script description and usage instructions.
# This script is used to calculate iHS from a VCF file.
# Usage: Rscript this_script.R input.vcf.gz output.csv

# Processing command line arguments
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]  # Path to the VCF file
output <- args[2]    # Path to the output file

# Reading the VCF file as a data table
vcf <- data.table::fread(vcf_file)

# Extracting an ID based on a specific position
id <- vcf %>%
  filter(POS == "500000") %>%
  select(ID)

# Loop to process chromosomes. Currently set for a single iteration.
for(i in 1:1) {
  # Generating the haplotype file name for each chromosome
  hap_file = paste("hap_chr_", i, ".cgu", sep = "")
  
  # Creating an internal representation of the haplotype data
  hh <- data2haplohh(hap_file = vcf_file,
                     polarize_vcf = FALSE,
                     chr.name = i,
                     vcf_reader = "data.table")
  
  # Performing the scan on a single chromosome to calculate iHH values
  scan <- scan_hh(hh)
  
  # Initializing or concatenating scan results for whole-genome analysis
  if (i == 1) {
    wgscan <- scan
  } else {
    wgscan <- rbind(wgscan, scan)
  }
}

# Calculating iHS from iHH values with frequency binning
wgscan.ihs <- ihh2ihs(wgscan, freqbin = 0.05)
ihs <- wgscan.ihs$ihs

# Preparing the output dataframe
ihs_output <- ihs %>%
  na.omit() %>%
  tibble::rownames_to_column("id") %>%
  filter(POSITION == 500000) %>%
  mutate(simulation = tools::file_path_sans_ext(vcf_file, compression = TRUE), vcf_ID = id) %>%
  select(simulation, vcf_ID, IHS)

# Creating the output file and appending the results
file.create(output)
write_tsv(ihs_output, output, col_names = TRUE, append = TRUE)


