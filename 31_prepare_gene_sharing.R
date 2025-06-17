# Author: dlu @ veelab
# Version: 2025-06-17

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/lisc/scratch/dome/luecking/projects/wwtp_linking")


# load blastp output ------------------------------------------------------

blast_out <- fread("intermediate/blastp/all_vs_all_blastp.tsv")
# blast_out <- fread("intermediate/blastp/all_vs_all_blastp_subset.tsv")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")


# add contig ids ----------------------------------------------------------

blast_out$q_contig_id <- str_remove(blast_out$qseqid, "\\_\\d*$")
blast_out$s_contig_id <- str_remove(blast_out$sseqid, "\\_\\d*$")


# filtering ---------------------------------------------------------------

blast_out <- blast_out %>% 
  filter(q_contig_id != s_contig_id) %>% 
  filter(pident >= 80, qcovhsp >= 80)



# summarize and save ------------------------------------------------------

summary_table <- blast_out %>% 
  group_by(q_contig_id, s_contig_id) %>% 
  summarise(genes_shared = n(), .groups = "drop") %>% 
  rename(from = q_contig_id, to = s_contig_id)


# count genes for each contig ---------------------------------------------

protein_file <- "intermediate/proteins/all_proteins.faa"
protein_accessions <- system2("grep", args = c("'>'", protein_file), stdout = TRUE)
protein_accessions <- as.data.table(protein_accessions)

# get to contig_id
protein_accessions$contig_id <- str_remove(str_remove(str_remove(protein_accessions$protein_accessions, "\\s.*$"), "\\_\\d*$"), ">")

protein_numbers <- as.data.table(table(protein_accessions$contig_id))


# calculate the gene sharing value ----------------------------------------

summary_table$genes_per_contig <- protein_numbers$N[match(summary_table$from, protein_numbers$V1)]
summary_table$gene_sharing <- summary_table$genes_shared / summary_table$genes_per_contig

fwrite(summary_table %>% select(from, to, gene_sharing), file = "intermediate/network/gene_sharing.csv")
