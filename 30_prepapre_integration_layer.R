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


# create two layers for each overhang type --------------------------------

blast_out8080$overhang_type <- if_else(
  stringr::str_ends(blast_out8080$qseqid, "m"), 
  true = "m", 
  false = "b"
)

for(type in unique(blast_out8080$overhang_type)){
  summary_table <- blast_out8080 %>%
    filter(overhang_type == type) %>% 
    group_by(mapped_to, sseqid) %>%
    summarise(number_of_reads = n(), .groups = 'drop') %>%
    rename(from = mapped_to, to = sseqid)
  
  summary_table$log10_reads <- log10(summary_table$number_of_reads + 1)
  names(summary_table) <- c("from", "to", "number_of_reads", paste0("integration_", type))
  fwrite(summary_table %>% select(from, to, paste0("integration_", type)), paste0("intermediate/network/integration_", type, ".csv"))
}



# explore -----------------------------------------------------------------

b_out <- blast_out %>% 
  filter(!str_ends(qseqid, "m")) %>% 
  filter(pident >= 70, qcov >= 5) %>% 
  mutate(direction_of_overhang = str_extract(qseqid, "\\d$"))

result_table <- b_out %>%
  group_by(mapped_to, sseqid) %>%
  summarise(
    number_of_5 = sum(direction_of_overhang == "5", na.rm = TRUE),
    number_of_3 = sum(direction_of_overhang == "3", na.rm = TRUE),
    .groups = 'drop' # This drops the grouping after summarising
  )

result_table_filtered <- result_table %>% 
  filter(number_of_5 > 0 & number_of_3 > 0)
