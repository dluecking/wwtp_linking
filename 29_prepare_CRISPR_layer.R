# Author: dlu @ veelab
# Version: 2025-06-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load data, clean ids ----------------------------------------------------

crispr_df <- fread("intermediate/CRISPR/map_spacers_to_targets/minced_results_spacers_filtered_lc_vs_target_db_all_hits.tsv")
crispr_df$from <- str_remove(crispr_df$`Spacer id`, "\\_CRISPR.*$")
crispr_df$to <- str_remove(crispr_df$`Target id`, "\\_polypolish$")


# filter self hits --------------------------------------------------------

crispr_df <- crispr_df %>% 
  filter(from != to)

summary_table <- crispr_df %>% 
  group_by(from, to) %>% 
  summarise(number_of_spacer_hits = n(), .groups = "drop")

summary_table$log10_number_of_hits <- log10(summary_table$number_of_spacer_hits + 1)


# save --------------------------------------------------------------------

names(summary_table) <- c("from", "to", "number_of_spacer_hits", "crispr") 
fwrite(summary_table, "intermediate/network/crispr.csv")

