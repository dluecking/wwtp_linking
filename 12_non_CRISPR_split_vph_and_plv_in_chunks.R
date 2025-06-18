# Author: dlu @ veelab
# Version: 2025-06-18

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# CONSTANTS
WINDOW <- 40
SLIDING <- 10

# load fasta --------------------------------------------------------------

fasta_vph <- lapply(list.files("intermediate/contigs/vph/", full.names = TRUE), read.fasta)
fasta_plv <- lapply(list.files("intermediate/contigs/plv/", full.names = TRUE), read.fasta)

fasta <- c(fasta_plv, fasta_vph)


#  write to chunks of WINDOW size -----------------------------------------

for(rec in fasta) {
  seq <- unlist(getSequence(rec, as.string = TRUE))
  starts <- seq(1, nchar(seq), by = SLIDING)
  chunks <- sapply(starts, function(i) substring(seq, i, min(i + (WINDOW-1), nchar(seq))))
  
  chunk_names <- paste0(getName(rec), "_chunk_", seq(1, length(chunks)))
  
  # Prepare all lines at once
  to_write <- paste0(">", chunk_names, "\n", chunks, collapse = "\n")
  
  # Write
  FILE <- paste0("intermediate/non_CRISPR/chunks/", getName(rec), "_chunks.fa")
  write(x = to_write, file = FILE)
}


