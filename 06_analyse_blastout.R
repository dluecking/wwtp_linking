# Author: dlu @ veelab
# Version: 2025-05-28

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load blast_out ----------------------------------------------------------

blast_out <- fread("intermediate/blast_results/query_vs_gv_contigs_megablast.tsv",
                   col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                 "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                                 "qlen", "slen"))

# Filter and process hits -------------------------------------------------

final_table_of_hits <- blast_out %>%
  as_tibble() %>% # Convert to tibble for consistent dplyr behavior
  mutate(
    pident = as.numeric(pident),
    qlen = as.numeric(qlen),
    length = as.numeric(length),
    qcov_hsp_perc = (length / qlen) * 100 # Calculate query coverage
  ) %>%
  filter(
    pident > 90, # Identity greater than 90%
    qcov_hsp_perc > 90, # Query coverage greater than 90%
    qseqid != sseqid # Exclude self-hits
  ) %>%
  select(
    query_contig_id = qseqid,
    gv_contig_id = sseqid,
    percent_identity = pident,
    query_coverage_percent = qcov_hsp_perc,
    evalue,
    bitscore,
    alignment_length = length,
    query_length = qlen,
    subject_length = slen
  )


# save to helperfile ------------------------------------------------------

fwrite(final_table_of_hits %>% select(query_contig_id), "helperfiles/contigs_to_remove_since_matching_GVs.txt", col.names = F)
