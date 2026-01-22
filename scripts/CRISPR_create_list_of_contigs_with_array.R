# Author: dlu @ veelab
# Version: 2025-06-10

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(ape)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load minced output ------------------------------------------------------

minced_gff <- read.gff("intermediate/minced/minced_results.gff") %>% 
  filter(type == "repeat_region")

nr_seqids <- minced_gff %>% 
  select(seqid) %>% 
  unique()

# Filter and write CRISPR candidate contigs -------------------------------

all_contigs <- read.fasta("intermediate/contigs/all_contigs.fna", seqtype = "DNA", as.string = TRUE)

crispr_contig_names_to_keep <- names(all_contigs)[names(all_contigs) %in% nr_seqids$seqid]

candidate_crispr_contigs <- all_contigs[crispr_contig_names_to_keep]

output_dir <- "intermediate/CRISPR"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_fasta_path <- file.path(output_dir, "candidate_contigs.fna")

write.fasta(
  sequences = getSequence(candidate_crispr_contigs),
  names = getName(candidate_crispr_contigs),
  file.out = output_fasta_path,
  nbchar = 60,
  open = "w"
)

message(paste("Filtered", length(candidate_crispr_contigs), "contigs with CRISPR arrays to:", output_fasta_path))
