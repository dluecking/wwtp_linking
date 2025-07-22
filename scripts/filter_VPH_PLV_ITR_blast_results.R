# Author: dlu @ veelab
# Version: 2025-07-18

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)


# deal with command line inputs -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript filter_VPH_PLV_ITR_blast_results.R <blast_results.in> <output.tsv>")
}

input_file  <- args[1]
output_file <- args[2]

if (!file.exists(input_file)) {
  stop(paste("Input file does not exist:", input_file))
}


# if debug
DEBUG <- TRUE
if(DEBUG){
  # set working directory
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  input_file <- "../intermediate/ITRs/ends_vs_ends.out"
  output_file <- "../intermediate/ITRs/vph_plv_ITR_info.tsv"
}


# load blast data ---------------------------------------------------------

blast_df <- fread(input_file)
names(blast_df) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                     "gapopen", "qstart", "qend", "sstart", "send",
                     "evalue", "bitscore", "qcovs")

blast_df$q_contig <- str_remove(blast_df$qseqid, "\\_[A-z]*$")
blast_df$s_contig <- str_remove(blast_df$sseqid, "\\_[A-z]*$")
blast_df <- blast_df %>% 
  filter(length >100) %>% # only hits longer than 100bp are valid
  filter(q_contig == s_contig & qseqid != sseqid) # only start vs end is valid but not self hits



# save data ---------------------------------------------------------------

fwrite(blast_df, file = output_file)
