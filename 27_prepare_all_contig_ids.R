# Author: dlu @ veelab
# Version: 2025-06-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load all contigs fasta --------------------------------------------------

fasta <- read.fasta("intermediate/contigs/all_contigs.fna")
names <- getName(fasta)
fwrite(names, "intermediate/network/all_contig_ids.txt")
