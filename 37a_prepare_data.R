# Author: dlu @ veelab
# Version: 2025-06-18
# Script: 37a_network_data_preparation.R
# Purpose: Load and prepare network data for all downstream analyses

# Packages ----------------------------------------------------------------
library(dplyr)
library(data.table)
library(stringr)
library(googlesheets4)
library(tidyr)

# Setup -------------------------------------------------------------------
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Create output directory if needed
dir.create("intermediate/network/network_analysis", recursive = TRUE, showWarnings = FALSE)


# Load contig information -------------------------------------------------
contig_df <- fread("intermediate/network/all_contig_ids.txt", header = F)
names(contig_df) <- c("contig_id")
contig_df$contig_id <- str_remove(contig_df$contig_id, ".polypolish$")
contig_df <- contig_df %>%
  mutate(
    type = case_when(
      str_ends(contig_id, "vph") ~ "vph",
      str_ends(contig_id, "lc")  ~ "lc",
      str_ends(contig_id, "plv") ~ "plv",
      TRUE ~ "gv"
    )
  )

# Add color scheme
contig_df <- contig_df %>% 
  mutate(
    color = case_when(
      str_ends(contig_id, "vph") ~ "goldenrod1",
      str_ends(contig_id, "lc")  ~ "seagreen",
      str_ends(contig_id, "plv") ~ "hotpink",
      TRUE ~ "steelblue"
    )
  )


# Load taxonomy information -----------------------------------------------
# GV info
sheet_url <- "https://docs.google.com/spreadsheets/d/1QLNiqSt0XOS4xVPAeZAppwVjjjPIKdEE6w6f2_Qm55c/edit?gid=1228834474#gid=1228834474"
GV_info <- read_sheet(sheet_url, sheet = "Final GVs overview") %>% 
  select(shortname, personal_assessment_order, public_ID) %>% 
  drop_na(shortname)

# VPH and PLV info
sheet_url <- "https://docs.google.com/spreadsheets/d/1CnqcfhOfS0rVBU6mvyCKrm50BmxJjVvif47f9vkQZgU/edit?gid=836958150#gid=836958150"
vph_info <- read_sheet(sheet_url, sheet = "Table S3") %>% 
  select(contig_ID, public_ID)
plv_info <- read_sheet(sheet_url, sheet = "Table S4") %>% 
  select(contig_ID, public_ID)
vph_plv_combined_info <- rbind(vph_info, plv_info)

# LC taxonomy info
lc_tax_info <- readRDS("intermediate/lc_tax/lc_tax_info_df.csv")
lc_tax_info_short <- rbindlist(lapply(list.files("intermediate/lc_tax/", pattern = ".*filtered.csv", full.names = TRUE), fread))


# Add taxonomy information to contig_df -----------------------------------
# First all lcs (with percentage tag)
contig_df$tax_info <- paste0(lc_tax_info_short$majority_organism[match(contig_df$contig_id, lc_tax_info_short$contig_id)],
                             " (", round(lc_tax_info_short$pct_reads_assigned_to_majority_taxon[match(contig_df$contig_id, lc_tax_info_short$contig_id)], 2), " %)")

# Then assign type-specific info
contig_df <- contig_df %>% 
  mutate(
    tax_info = case_when(
      str_ends(contig_id, "vph") ~ contig_id,
      str_ends(contig_id, "plv") ~ contig_id,
      str_ends(contig_id, "lc") ~ tax_info,
      TRUE ~ "GV - to be added"
    )
  )

# Add GV taxonomy
for(i in 1:nrow(contig_df)){
  if(contig_df$type[i] == "gv")
    if(length((GV_info$personal_assessment_order[GV_info == contig_df$contig_id[i]])) > 0){
      contig_df$tax_info[i] <- GV_info$personal_assessment_order[GV_info == contig_df$contig_id[i]]
    } else {
      contig_df$tax_info[i] <- "unknown"
    }
}


# Load CRISPR data --------------------------------------------------------
crispr_df <- fread("intermediate/CRISPR/cassette/HMM2019_cassettes.csv") %>% 
  filter(bitscore >= 100)
crispr_df$contig_id <- str_remove(crispr_df$V1, "\\_\\d+\\_ID.*$")


# Load edge lists ---------------------------------------------------------
edgelist_crispr <- fread("intermediate/network/crispr.csv")
edgelist_integration_b <- fread("intermediate/network/integration_b.csv")
edgelist_integration_m <- fread("intermediate/network/integration_m.csv")
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing.csv")
edgelist_occurance_ill <- fread("intermediate/network/occurance_ill_strat.csv") 
edgelist_occurance_ont <- fread("intermediate/network/occurance_ont_strat.csv")
edgelist_non_crispr <- fread("intermediate/network/non_CRISPR.csv")


# load membership list ----------------------------------------------------

cat("\n=== Loading LC cluster mapping ===\n")
# membership_df <- fread("intermediate/mash/lc_contigs_only_cluster_membership.csv")
membership_df <- fread("intermediate/mash/lc_contigs_only_cluster_membership_005.csv")

# Process occurrence edges (Illumina and ONT) ------------------------------
# Add edge IDs and clean Illumina df
edgelist_occurance_ill <- edgelist_occurance_ill %>%
  mutate(
    node1 = pmin(from, to),
    node2 = pmax(from, to),
    edge_id = paste0(node1, "--", node2)
  ) %>%
  select(edge_id, node1, node2, 
         spearman_ill = spearman_ont, # this is a mistake in another script
         correlation_ill = correlation) %>%
  distinct(edge_id, .keep_all = TRUE)

# Same for ONT
edgelist_occurance_ont <- edgelist_occurance_ont %>%
  mutate(
    node1 = pmin(from, to),
    node2 = pmax(from, to),
    edge_id = paste0(node1, "--", node2)
  ) %>%
  select(edge_id, node1, node2, 
         spearman_ont = spearman_ont, 
         correlation_ont = correlation) %>%
  distinct(edge_id, .keep_all = TRUE)

# Inner join to get only shared edges
occurance_edges <- edgelist_occurance_ill %>%
  inner_join(edgelist_occurance_ont, by = c("edge_id", "node1", "node2")) %>%
  mutate(
    agree = (correlation_ill == correlation_ont)
  ) %>%
  select(edge_id, node1, node2, 
         spearman_ill, correlation_ill, 
         spearman_ont, correlation_ont, 
         agree)

# Remove disagreements
occurance_edges <- occurance_edges %>% filter(agree == "TRUE")

# Create value string
occurance_edges$value <- paste0(round(occurance_edges$spearman_ill, 2), "_ill--", round(occurance_edges$spearman_ont, 2), "_ont")

# Rename columns
names(occurance_edges) <- c("edge_id", "from", "to", "spearman_ill", "correlation_ill", "spearman_ont", "correlation_ont", "agree", "value") 


# Create combined edge dataframe ------------------------------------------

cat("\n=== Creating combined edge dataframe with weights ===\n")

# Define weight constants (adjust these to tune clustering behavior)
WEIGHT_CRISPR <- 1000
WEIGHT_GENE_SHARING <- 3
WEIGHT_INTEGRATION_BOUNDARY <- 1000
WEIGHT_INTEGRATION_MIDDLE <- 750
WEIGHT_NON_CRISPR <- 1000

# Co-occurrence weights by correlation strength
WEIGHT_COOCCUR_VERY_STRONG_POS <- 3   # Spearman >= 0.80, positive
WEIGHT_COOCCUR_STRONG_POS <- 2        # Spearman >= 0.70, positive
WEIGHT_COOCCUR_MODERATE_POS <- 1      # Spearman >= 0.60, positive

WEIGHT_COOCCUR_VERY_STRONG_NEG <- 3   # Spearman >= 0.80, negative
WEIGHT_COOCCUR_STRONG_NEG <- 2        # Spearman >= 0.70, negative
WEIGHT_COOCCUR_MODERATE_NEG <- 1      # Spearman >= 0.60, negative

# Spearman thresholds for co-occurrence weight tiers
SPEARMAN_VERY_STRONG <- 0.80
SPEARMAN_STRONG <- 0.70

# Calculate average Spearman for co-occurrence edges
occurance_edges <- occurance_edges %>%
  mutate(
    spearman_avg = (abs(spearman_ill) + abs(spearman_ont)) / 2
  )

big_connection_df <- rbind(
  # High-confidence mechanistic evidence
  edgelist_crispr %>% 
    select(from, to, value = "crispr") %>% 
    mutate(type = "crispr", weight = WEIGHT_CRISPR),
  
  edgelist_gene_sharing %>% 
    select(from, to, value = "gene_sharing") %>% 
    mutate(type = "gene_sharing", weight = WEIGHT_GENE_SHARING),
  
  edgelist_integration_b %>% 
    select(from, to, value = "integration_b") %>% 
    mutate(type = "integration_boundary", weight = WEIGHT_INTEGRATION_BOUNDARY),
  
  edgelist_integration_m %>% 
    select(from, to, value = "integration_m") %>% 
    mutate(type = "integration_middle", weight = WEIGHT_INTEGRATION_MIDDLE),
  
  # Medium-confidence evidence
  edgelist_non_crispr %>% 
    select(from, to, value = "non_CRISPR") %>% 
    mutate(type = "non_crispr", weight = WEIGHT_NON_CRISPR),
  
  # Positive co-occurrence (weighted by correlation strength)
  occurance_edges %>%
    filter(correlation_ont == "positive") %>% 
    mutate(
      type = "occurance_positive",
      weight = case_when(
        spearman_avg >= SPEARMAN_VERY_STRONG ~ WEIGHT_COOCCUR_VERY_STRONG_POS,
        spearman_avg >= SPEARMAN_STRONG ~ WEIGHT_COOCCUR_STRONG_POS,
        TRUE ~ WEIGHT_COOCCUR_MODERATE_POS
      )
    ) %>%
    select(from, to, value, type, weight),
  
  # Negative co-occurrence (lower weights)
  occurance_edges %>%
    filter(correlation_ont == "negative") %>% 
    mutate(
      type = "occurance_negative",
      weight = case_when(
        spearman_avg >= SPEARMAN_VERY_STRONG ~ WEIGHT_COOCCUR_VERY_STRONG_NEG,
        spearman_avg >= SPEARMAN_STRONG ~ WEIGHT_COOCCUR_STRONG_NEG,
        TRUE ~ WEIGHT_COOCCUR_MODERATE_NEG
      )
    ) %>%
    select(from, to, value, type, weight)
)

cat("Edge weights configured:\n")
cat("  Mechanistic evidence (CRISPR/integration/genes):", WEIGHT_CRISPR, "\n")
cat("  Non-CRISPR evidence:", WEIGHT_NON_CRISPR, "\n")
cat("  Co-occurrence (very strong ≥", SPEARMAN_VERY_STRONG, "):", 
    WEIGHT_COOCCUR_VERY_STRONG_POS, "(pos),", WEIGHT_COOCCUR_VERY_STRONG_NEG, "(neg)\n")
cat("  Co-occurrence (strong ≥", SPEARMAN_STRONG, "):", 
    WEIGHT_COOCCUR_STRONG_POS, "(pos),", WEIGHT_COOCCUR_STRONG_NEG, "(neg)\n")
cat("  Co-occurrence (moderate):", 
    WEIGHT_COOCCUR_MODERATE_POS, "(pos),", WEIGHT_COOCCUR_MODERATE_NEG, "(neg)\n\n")


# add membership information ----------------------------------------------
# remove polypolish from GV name
membership_df$contig_id <- str_remove(membership_df$contig_id, "\\spolypolish")
membership_df$cluster_id <- str_remove(membership_df$cluster_id, "\\spolypolish")

# replace from to with the cluster name (which is the same for GVs, PLVs, VPHs, but different for LCs)
big_connection_df$from <- membership_df$cluster_id[match(big_connection_df$from, membership_df$contig_id)]
big_connection_df$to <- membership_df$cluster_id[match(big_connection_df$to, membership_df$contig_id)]


# Create unique edge IDs and filter duplicates ----------------------------
big_connection_df_filtered <- big_connection_df %>%
  mutate(
    # For undirected edges, create sorted edge_id
    edge_id_undirected = paste0(pmin(from, to), "--", pmax(from, to)),
    # For directed edges (CRISPR, non_crispr), keep original order
    edge_id = if_else(
      type %in% c("crispr", "non_crispr"),
      paste0(from, "--", to),  # Keep direction
      edge_id_undirected       # Sort alphabetically
    ),
    # Mark which edges are directed
    is_directed = type %in% c("crispr", "non_crispr")
  ) %>%
  select(-edge_id_undirected)

# Keep only one row per edge_id + type
big_connection_df_filtered <- big_connection_df_filtered %>% 
  dplyr::distinct(edge_id, type, .keep_all = TRUE)

# remove "self connections" due to clustering:
big_connection_df_filtered <- big_connection_df_filtered %>% filter(from != to)

# quick testing
# for(cluster in a %>% arrange(desc(N)) %>% top_n(10) %>% pull(V1)){
#   cat(cluster, "\n")
#   cat("occurs this many times UNFILTERED and uncollapsed:\n")
#   s <- big_connection_df %>% 
#     filter(if_any(everything(), ~ grepl(cluster, .))) %>% 
#     nrow()
#   print(s)
#   cat("occurs this many times FILTERED and collapsed:\n")
#   s <- big_connection_df_filtered %>% 
#     filter(if_any(everything(), ~ grepl(cluster, .))) %>% 
#     nrow()
#   print(s)
# }



# Add node type information to edges --------------------------------------
big_connection_df_filtered <- big_connection_df_filtered %>% 
  mutate(from_type = case_when(
    str_starts(from, "LC_CLUSTER")  ~ "lc",
    str_ends(from, "vph") ~ "vph",
    str_ends(from, "plv") ~ "plv",
    TRUE ~ "gv"
  )) %>% 
  mutate(to_type = case_when(
    str_starts(to, "LC_CLUSTER")  ~ "lc",
    str_ends(to, "vph") ~ "vph",
    str_ends(to, "plv") ~ "plv",
    TRUE ~ "gv"
  ))


# Save prepared data ------------------------------------------------------
# Save as CSV where possible (main dataframes)
fwrite(big_connection_df_filtered, "intermediate/network/network_analysis/big_connection_df_filtered.csv")
fwrite(contig_df, "intermediate/network/network_analysis/contig_df.csv")
fwrite(GV_info, "intermediate/network/network_analysis/GV_info.csv")
fwrite(vph_plv_combined_info, "intermediate/network/network_analysis/vph_plv_combined_info.csv")
fwrite(crispr_df, "intermediate/network/network_analysis/crispr_df.csv")
fwrite(lc_tax_info_short, "intermediate/network/network_analysis/lc_tax_info_short.csv")

cat("\n=== Data preparation complete ===\n")
cat("Files saved to: intermediate/network/network_analysis/\n")
cat("  - big_connection_df_filtered.csv (", nrow(big_connection_df_filtered), " edges)\n")
cat("  - contig_df.csv (", nrow(contig_df), " nodes)\n")
cat("  - GV_info.csv\n")
cat("  - vph_plv_combined_info.csv\n")
cat("  - crispr_df.csv\n")
cat("  - lc_tax_info_short.csv\n\n")

# Weight summary ----------------------------------------------------------
cat("\n=== Edge weight summary ===\n")
weight_summary <- big_connection_df_filtered %>%
  group_by(type) %>%
  summarize(
    n_edges = n(),
    mean_weight = round(mean(weight), 2),
    median_weight = median(weight),
    min_weight = min(weight),
    max_weight = max(weight)
  ) %>%
  arrange(desc(mean_weight))

print(weight_summary)

# Summary statistics
cat("=== Network Summary ===\n")
cat("Total edges:", nrow(big_connection_df_filtered), "\n")
cat("Total nodes:", nrow(contig_df), "\n")
cat("Edge types:\n")
print(table(big_connection_df_filtered$type))
cat("\nNode types:\n")
print(table(contig_df$type))
