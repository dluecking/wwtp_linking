# Author: dlu @ veelab
# Version: 2025-05-28

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load fasta --------------------------------------------------------------

gene_names <- getName(read.fasta("intermediate/proteins/all_proteins.faa"))

gene_to_genome_df <- data.table(protein_id = gene_names,
                                contig_id = str_remove(gene_names, "\\_\\d*$"),
                                keyword = "")
fwrite(gene_to_genome_df, "intermediate/proteins/all_proteins_gene_to_genome.csv", sep = ",")