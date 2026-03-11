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




# exploratory PAM motifs --------------------------------------------------

crispr_df <- crispr_df %>%
  rowwise() %>%
  mutate(pam_motif_present = {
    # PAM patterns to search
    pam_patterns <- c("NGG" = "[ATCG]GG", "NAG" = "[ATCG]AG", "NNGG" = "[ATCG]{2}GG")
    
    # Search region: last 5bp of upstream + first 5bp of downstream
    up_region <- substr(upstream, nchar(upstream) - 4, nchar(upstream))
    down_region <- substr(downstream, 1, 5)
    search_region <- paste0(up_region, down_region)
    
    # Find first matching PAM
    pam_found <- "none"
    for (pam_name in names(pam_patterns)) {
      if (str_detect(search_region, pam_patterns[pam_name])) {
        # Extract the actual motif
        match <- str_extract(search_region, pam_patterns[pam_name])
        pam_found <- match[1]
        break
      }
    }
    pam_found
  }) %>%
  ungroup()


# Filter for trusted hits (≤2 mismatches)
trusted_hits <- crispr_df %>%
  filter(`N mismatches` <= 3 & pam_motif_present != "none")

trusted_hits <- trusted_hits %>%
  mutate(
    source_is_NCV = !grepl("_lc", from),
    target_is_NCV = !grepl("_lc", to),
    interaction_type = case_when(
      source_is_NCV & target_is_NCV ~ "NCV -> NCV",
      source_is_NCV & !target_is_NCV ~ "NCV -> Microbial",
      !source_is_NCV & target_is_NCV ~ "Microbial -> NCV",
      TRUE ~ "Microbial -> Microbial"
    )
  )

# Detailed summary: from -> to with counts
detailed_summary <- trusted_hits %>%
  group_by(from, to, interaction_type) %>%
  summarise(
    n_matches = n(),
    .groups = "drop"
  ) %>%
  arrange(interaction_type, desc(n_matches))

# Summary by interaction type
type_summary <- detailed_summary %>%
  group_by(interaction_type) %>%
  summarise(
    n_connections = n(),
    total_matches = sum(n_matches),
    mean_matches = round(mean(n_matches), 1),
    median_matches = median(n_matches)
  )
