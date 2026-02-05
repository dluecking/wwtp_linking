# Author: dlu @ veelab
# Version: 2025-06-18
# Script: 37f_exploratory_network_stuff.R
# Purpose: Suggestions for additional analyses and virus-host prediction exploration

# Packages ----------------------------------------------------------------
library(dplyr)
library(data.table)
library(ggplot2)
library(igraph)

# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")

cat("\n")
cat("=====================================================\n")
cat("  SUGGESTED ANALYSES & VIRUS-HOST PREDICTIONS\n")
cat("=====================================================\n\n")


# =========================================================================
# VIRUS-HOST PREDICTION ANALYSIS
# =========================================================================

cat("=== VIRUS-HOST PREDICTION ANALYSIS ===\n\n")

# Count virus-LC connections by edge type
virus_host_candidates <- big_connection_df_filtered %>%
  filter(
    (from_type %in% c("vph", "plv", "gv") & to_type == "lc") |
      (to_type %in% c("vph", "plv", "gv") & from_type == "lc")
  ) %>%
  mutate(
    virus = ifelse(from_type != "lc", from, to),
    host = ifelse(from_type == "lc", from, to),
    virus_type = ifelse(from_type != "lc", from_type, to_type)
  ) %>%
  group_by(virus, host, virus_type) %>%
  summarise(
    n_edge_types = n(),
    edge_types = paste(type, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_edge_types))

cat("Virus-host pairs with multiple evidence types:\n")
print(head(virus_host_candidates, 20))

cat("\nSummary by virus type:\n")
print(virus_host_candidates %>% 
        group_by(virus_type) %>% 
        summarise(
          n_pairs = n(),
          pairs_multi_evidence = sum(n_edge_types > 1),
          max_evidence_types = max(n_edge_types)
        ))

# Save for further analysis
fwrite(virus_host_candidates, "testing/virus_host_candidates.csv")
cat("\nSaved to: testing/virus_host_candidates.csv\n")


# =========================================================================
# CO-OCCURRENCE NETWORK SUMMARY
# =========================================================================

cat("\n\n=== CO-OCCURRENCE NETWORK SUMMARY ===\n\n")

occurance_summary <- big_connection_df_filtered %>%
  filter(type %in% c("occurance_positive", "occurance_negative")) %>%
  mutate(
    interaction = paste0(
      pmin(from_type, to_type), 
      " <-> ", 
      pmax(from_type, to_type)
    )
  ) %>%
  group_by(type, interaction) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(type, desc(n))

cat("Co-occurrence interactions by node type:\n")
print(occurance_summary)

cat("\nSee virus_microbe_correlations.md for biological interpretation!\n")


# =========================================================================
# EDGE OVERLAP ANALYSIS
# =========================================================================

cat("\n\n=== EDGE OVERLAP ANALYSIS ===\n\n")

edge_overlap <- big_connection_df_filtered %>%
  mutate(edge_pair = paste0(pmin(from, to), "--", pmax(from, to))) %>%
  group_by(edge_pair) %>%
  summarise(
    n_edge_types = n(),
    edge_types = paste(type, collapse = ", ")
  ) %>%
  arrange(desc(n_edge_types))

cat("Edges with multiple types:\n")
print(table(edge_overlap$n_edge_types))

cat("\nTop edges with most types:\n")
print(head(edge_overlap %>% filter(n_edge_types > 1), 10))


# =========================================================================
# TRIANGLE MOTIF ANALYSIS
# =========================================================================

cat("\n\n=== NETWORK MOTIF ANALYSIS ===\n\n")

g <- graph_from_data_frame(big_connection_df_filtered, directed = FALSE)
triangles_count <- sum(count_triangles(g)) / 3  # Each triangle counted 3 times
cat("Number of triangles in full network:", triangles_count, "\n")


# =========================================================================
# K-CORE DECOMPOSITION
# =========================================================================

cat("\n\n=== CORE-PERIPHERY STRUCTURE ===\n\n")

coreness_vals <- coreness(g)
cat("K-core distribution:\n")
print(table(coreness_vals))

cat("\nHighest k-core value:", max(coreness_vals), "\n")


# =========================================================================
# ANALYSIS SUGGESTIONS
# =========================================================================

cat("\n\n")
cat("==========================================================\n")
cat("PRIORITY RECOMMENDATIONS:\n")
cat("==========================================================\n\n")

cat("HIGH PRIORITY:\n")
cat("--------------\n")
cat("1. Virus-host prediction validation\n")
cat("   - Use virus_host_candidates.csv\n")
cat("   - Score predictions by number of edge types\n")
cat("   - Validate against known host relationships\n")
cat("   - Compare CRISPR vs integration vs gene sharing\n\n")

cat("2. Co-occurrence network interpretation\n")
cat("   - Positive vs negative correlations\n")
cat("   - See virus_microbe_correlations.md for predictions\n")
cat("   - Test if patterns match biological expectations\n\n")

cat("3. Multi-layer network analysis\n")
cat("   - How do different edge types contribute?\n")
cat("   - Edge overlap between layers\n")
cat("   - Layer-specific centrality\n\n")

cat("4. Weighted network analysis\n")
cat("   - Gene sharing edges have similarity scores\n")
cat("   - Weighted centrality metrics\n")
cat("   - Minimum spanning tree\n\n")

cat("MEDIUM PRIORITY:\n")
cat("----------------\n")
cat("5. Network motif analysis\n")
cat("   - 3-node motifs (triads)\n")
cat("   - Over/under-represented patterns\n")
cat("   - Biological interpretation of triangles\n\n")

cat("6. Link prediction\n")
cat("   - Common neighbors score\n")
cat("   - Adamic-Adar index\n")
cat("   - Predict missing virus-host relationships\n\n")

cat("7. Functional enrichment\n")
cat("   - Test if highly connected LCs belong to specific taxa\n")
cat("   - Test if GV-connected LCs have specific profiles\n\n")

cat("LOWER PRIORITY:\n")
cat("---------------\n")
cat("8. Backbone extraction (if network too large)\n")
cat("9. Overlapping communities (link communities)\n")
cat("10. Modularity optimization across resolutions\n")
cat("11. Differential network analysis (if multiple conditions)\n")
cat("12. Permutation tests for significance\n\n")

cat("==========================================================\n\n")


# =========================================================================
# SPECIFIC CODE TEMPLATES
# =========================================================================

cat("CODE TEMPLATES:\n")
cat("---------------\n\n")

cat("# Example: Filter high-confidence virus-host predictions\n")
cat("high_confidence <- virus_host_candidates %>%\n")
cat("  filter(n_edge_types >= 2) %>%\n")
cat("  filter(grepl('crispr|integration', edge_types))\n\n")

cat("# Example: Build positive co-occurrence network only\n")
cat("positive_net <- big_connection_df_filtered %>%\n")
cat("  filter(type == 'occurance_positive')\n")
cat("g_pos <- graph_from_data_frame(positive_net, directed = FALSE)\n\n")

cat("# Example: Calculate weighted centrality for gene sharing\n")
cat("gene_sharing_net <- big_connection_df_filtered %>%\n")
cat("  filter(type == 'gene_sharing')\n")
cat("g_gs <- graph_from_data_frame(gene_sharing_net, directed = FALSE)\n")
cat("E(g_gs)$weight <- as.numeric(gene_sharing_net$value)\n")
cat("weighted_betw <- betweenness(g_gs, weights = E(g_gs)$weight)\n\n")

cat("==========================================================\n")
cat("See individual sections above for implementation ideas!\n")
cat("==========================================================\n\n")

# Clean up
rm(g)
