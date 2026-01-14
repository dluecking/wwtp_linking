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
final_table_of_hits <- blast_out %>%
  as_tibble() %>%
  mutate(
    pident = as.numeric(pident),
    qlen = as.numeric(qlen),
    length = as.numeric(length)
  ) %>%
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
    total_query_coverage = (total_alignment_length / query_length) * 100,
    total_subject_coverage = (total_alignment_length / subject_length) * 100
  ) %>%
  # Now you can filter effectively
  filter(
    avg_identity > 90,
    total_query_coverage > 90 | total_subject_coverage > 90
  ) %>% 
  # now lets only keep same sample hits
  mutate(
    q_sample = str_remove(qseqid, "\\_.*$"),
    s_sample = str_remove(sseqid, "\\_.*$")
  ) %>% 
  filter(q_sample == s_sample)



# save to helperfile ------------------------------------------------------

fwrite(final_table_of_hits %>% select(qseqid) %>% unique(), "helperfiles/contigs_to_remove_since_matching_GVs.txt", col.names = F)
