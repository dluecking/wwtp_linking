# Author: dlu @ veelab
# Version: 2026-07-23

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(cowplot)

print("running script in review mode!")
# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")

# what fraction of genes needs to be shared so I keep it?
CUTOFF <- 0.05 # this is 5% right now

# load blastp output ------------------------------------------------------

blast_out <- fread("intermediate/blastp/all_vs_all_blastp.tsv")
# blast_out <- fread("../intermediate/blastp/all_vs_all_blastp_subset.tsv")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")


# add contig ids ----------------------------------------------------------

blast_out$q_contig_id <- str_remove(blast_out$qseqid, "\\_\\d*$")
blast_out$s_contig_id <- str_remove(blast_out$sseqid, "\\_\\d*$")


# filtering ---------------------------------------------------------------

blast_out <- blast_out %>% 
  filter(q_contig_id != s_contig_id) %>% 
  filter(pident >= 80, qcovhsp >= 80) %>% 
  group_by(qseqid, sseqid) %>% 
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>% 
  ungroup()


# REVIEWED NEW TABLE ------------------------------------------------------

# instead of getting to: contig_id_A, contig_id_B, percentage_of_genes_shared we aim for:
# contig_id_A, contig_id_B, list_of_pidents_of_genes_shared_above8080, total number of genes the smaller contig has

# step 1: get shared gene pidents into a list:
summary_table <- blast_out %>% 
  group_by(q_contig_id, s_contig_id) %>% 
  summarise(
    genes_shared_count = n(),
    pident_values = list(pident),
    .groups = "drop"
  ) %>% 
  rename(from = q_contig_id, to = s_contig_id)


# step 2: count genes for each contig

protein_file <- "intermediate/proteins/all_proteins.faa"
protein_accessions <- system2("grep", args = c("'>'", protein_file), stdout = TRUE)
protein_accessions <- as.data.table(protein_accessions)

# get to contig_id
protein_accessions$contig_id <- str_remove(str_remove(str_remove(protein_accessions$protein_accessions, "\\s.*$"), "\\_\\d*$"), ">")
protein_numbers <- as.data.table(table(protein_accessions$contig_id))

# step 3: add gene count per genome
summary_table$genes_per_contig_FROM <- protein_numbers$N[match(summary_table$from, protein_numbers$V1)]
summary_table$genes_per_contig_TO <- protein_numbers$N[match(summary_table$to, protein_numbers$V1)]


# step 4: calculate new weights:
summary_table <- summary_table %>% 
  mutate(
    # how many genes does the smaller contig of the connection hold?
    min_contig_genes = pmin(genes_per_contig_FROM, genes_per_contig_TO),
    
    # percentage of genes shared of this smaller contig
    prop_shared = genes_shared_count / min_contig_genes
    ) %>% 
  
  # remove hits where we have less than 5% genes shared
  filter(prop_shared >= CUTOFF) %>% 
  
  mutate(
    # double proof that we dont share more than 100%
    prop_shared = pmin(prop_shared, 1.0),
    
    # calc mean of the shared genes
    mean_pident_frac = sapply(pident_values, mean) / 100,
    
    # cal the raw score multiplying the number of genes with the mean pident
    raw_score = prop_shared * mean_pident_frac,
  )

fwrite(summary_table, file = "intermediate/network/REVIEW_gene_sharing.csv")