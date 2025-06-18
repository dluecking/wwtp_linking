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


# load data ---------------------------------------------------------------

blast_files <- list.files("intermediate/non_CRISPR/blastn_out/", full.names = TRUE)
blast_files <- blast_files[file.info(blast_files)$size > 0]

blast_out <- rbindlist(lapply(blast_files, fread), use.names = TRUE, fill = TRUE)


blast_out$to <- str_remove(blast_out$V2, "\\_chunk.*$")



# prepare outtable --------------------------------------------------------

summary_table <- blast_out %>% 
  group_by(V1, to) %>% 
  summarise(connections = n(), .groups = "drop")

names(summary_table) <- c("from", "to", "non_CRISPR")

fwrite(summary_table, "intermediate/network/non_CRISPR.csv")
