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

# this one is better, but not implemented yet:
final_table_of_hits <- blast_out_NEW %>%
  as_tibble() %>%
  mutate(
    pident = as.numeric(pident),
    qlen = as.numeric(qlen),
    length = as.numeric(length)
  ) %>%
  # GROUP BY QUERY: This is the critical missing step
  group_by(qseqid, sseqid) %>%
  summarise(
    total_alignment_length = sum(length),
    query_length = first(qlen),
    subject_length = first(slen),
    # We take the mean or max identity of the segments
    avg_identity = mean(pident),
    .groups = "drop"
  ) %>%
  mutate(
    total_query_coverage = (total_alignment_length / query_length) * 100
  ) %>%
  # Now you can filter effectively
  filter(
    avg_identity > 90,
    total_query_coverage > 90
  )

# save to helperfile ------------------------------------------------------

fwrite(final_table_of_hits %>% select(qseqid), "helperfiles/contigs_to_remove_since_matching_GVs.txt", col.names = F)
