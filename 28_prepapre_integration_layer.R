# Author: dlu @ veelab
# Version: 2025-06-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load data and retrieve mapped_to ----------------------------------------

blast_out <- fread("intermediate/integration/blast_results/all_overhangs_blastn.tsv")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                      "qlen", "slen")

blast_out$read <- str_remove(blast_out$qseqid, "\\_to\\_.*$")
blast_out$mapped_to <- str_remove(str_remove(blast_out$qseqid, "^.*\\_to\\_"), "\\_overhang.*$")

# filter out self hits
blast_out <- blast_out %>% 
  filter(mapped_to != sseqid)

blast_out$qcov <- blast_out$length/blast_out$qlen * 100

blast_out8080 <- blast_out %>% filter(pident >= 80, qcov >= 80)

summary_table <- blast_out8080 %>%
  group_by(mapped_to, sseqid) %>%
  summarise(number_of_reads = n(), .groups = 'drop') %>%
  rename(from = mapped_to, to = sseqid)

summary_table$log10_reads <- log10(summary_table$number_of_reads + 1)

names(summary_table) <- c("from", "to", "number_of_reads", "integration")
fwrite(summary_table %>% select(from, to, integration), "intermediate/network/integration.csv")
