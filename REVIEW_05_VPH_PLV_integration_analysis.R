# Author: dlu @ veelab
# Version: 2026-07-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load blast out ----------------------------------------------------------

blast_out <- fread("intermediate/REVIEW_vph_plv_integration/vph_plv_vs_lc_megablast.tsv")
