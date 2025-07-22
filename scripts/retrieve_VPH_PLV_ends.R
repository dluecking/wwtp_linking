# Author: dlu @ veelab
# Version: 2025-07-18

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)


# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# deal with input ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript scripts/retrieve_VPH_PLV_ends.R <output_file>")
}

OUTPUT_FILE <- args[1]


# load fasta --------------------------------------------------------------
VPH_DIR <- "intermediate/contigs/vph"
PLV_DIR <- "intermediate/contigs/plv"
vph_fasta <- lapply(list.files(VPH_DIR, full.names = TRUE), read.fasta)
plv_fasta <- lapply(list.files(PLV_DIR, full.names = TRUE), read.fasta)

fasta <- c(vph_fasta, plv_fasta)


WINDOW <- 1000

for(record in fasta){
  seq <- unlist(getSequence(record, as.string = TRUE))
  acc <- getName(record)

  start <- substring(seq, 1, WINDOW)
  end <- substring(seq, nchar(seq)-WINDOW+1, nchar(seq))

  to_write <- paste0(">", acc, "_start\n",
                     start, "\n",
                     ">", acc, "_end\n",
                     end)
  write(to_write, file = OUTPUT_FILE, append = TRUE)

}