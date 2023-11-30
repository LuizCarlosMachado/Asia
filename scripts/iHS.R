library(rehh)
library(dplyr)
library(readr)


# Script to calculate iHS
# usage: Rscript input.vcf.gz derived.allele_position output.csv
# DATABASE
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output <- args[2]


vcf <- data.table::fread( vcf_file)
id <- vcf %>%
  filter(POS == "500000") %>%
  select(ID) 


for(i in 1:1) {
  # haplotype file name for each chromosome
  hap_file = paste("hap_chr_", i, ".cgu", sep = "")
  # create internal representation
  hh <- data2haplohh(hap_file = vcf_file,
                     polarize_vcf = FALSE,
                     chr.name = i,
                     vcf_reader = "data.table")
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    wgscan <- scan
  } else {
    wgscan <- rbind(wgscan, scan)
  }
}

wgscan.ihs <- ihh2ihs(wgscan, freqbin = 0.05)
ihs <- wgscan.ihs$ihs

ihs_output <- ihs %>%
  na.omit() %>%
  tibble::rownames_to_column("id") %>%
  filter(POSITION == 500000) %>%
  mutate(simulation = tools::file_path_sans_ext(vcf_file, compression = T), vcf_ID = id) %>%
  select(simulation, vcf_ID, IHS)

file.create(output)
write_tsv(ihs_output, output, col_names = TRUE, append = TRUE)

