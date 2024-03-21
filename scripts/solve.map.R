# Load the necessary libraries
library(readr)  # For reading and writing tabular data
library(purrr)  # For functional programming tools, though it's not directly used in the script

# Parse command line arguments to get the input and output file paths
args <- commandArgs(trailingOnly = TRUE)
map_file <- args[1]  # The path to the input data file
output <- args[2]    # The path to the output file

# Read the input file and modify it
map <- data.table::fread(map_file) %>%
  dplyr::mutate(V3 = V4/10^6)  # Modify the third column of the data by dividing the fourth column by 10^6

# Create the output file if it does not exist
file.create(output)

# Write the modified data to the output file, without column names and allowing for appending
write_tsv(map, output, col_names = FALSE, append = TRUE)