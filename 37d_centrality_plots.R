# Author: dlu @ veelab
# Version: 2025-06-18
# Script: 37d_centrality_analysis.R
# Purpose: Calculate and visualize network centrality metrics (degree, betweenness)
#          with special focus on LC contigs connected to GVs

# Packages ----------------------------------------------------------------
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(igraph)
library(cowplot)
library(ggsignif)
library(tibble)
library(doParallel)
library(foreach)
library(patchwork)

# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# Parallel setup
cores_to_use <- 14
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# Create output directories
dir.create("final/centrality_plots", recursive = TRUE, showWarnings = FALSE)


# Load prepared data ------------------------------------------------------
big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")
membership_df <- fread("intermediate/mash/lc_contigs_only_cluster_membership_005.csv")
# Remove polypolish suffix from membership
membership_df$contig_id <- str_remove(membership_df$contig_id, "\\spolypolish")
membership_df$cluster_id <- str_remove(membership_df$cluster_id, "\\spolypolish")

# Create nodes_df for graph construction
nodes_df <- contig_df %>%
  left_join(membership_df %>% select(contig_id, cluster_id, cluster_size), 
            by = "contig_id") %>%
  # Use cluster_id where available, otherwise keep original
  mutate(
    node_id = ifelse(!is.na(cluster_id), cluster_id, contig_id),
    original_contig_id = contig_id  # Keep original
  ) %>%
  # CRITICAL: Keep only ONE row per node_id (first occurrence)
  distinct(node_id, .keep_all = TRUE) %>%
  # Rename for igraph compatibility
  select(-contig_id) %>%
  rename(contig_id = node_id)

# Update type for LC clusters
nodes_df <- nodes_df %>%
  mutate(
    type = case_when(
      grepl("LC_CLUSTER", contig_id) ~ "lc",
      str_ends(contig_id, "vph") ~ "vph",
      str_ends(original_contig_id, "lc")  ~ "lc",
      str_ends(contig_id, "plv") ~ "plv",
      TRUE ~ "gv"
    )
  )

# make sure contig_id is early:
nodes_df <- nodes_df %>%
  select(contig_id, everything())


cat("Loaded data:\n")
cat("  - Edges:", nrow(big_connection_df_filtered), "\n")
cat("  - Nodes:", nrow(contig_df), "\n\n")


# Identify LC contigs connected to GVs ------------------------------------
cat("=== Identifying LC contigs connected to GVs ===\n")

lc_gv_connected <- big_connection_df_filtered %>%
  filter(from_type == "gv" | to_type == "gv") %>%   # Only GV edges
  mutate(lc = case_when(
    str_starts(from, "LC_CLUSTER") ~ from,
    str_starts(to, "LC_CLUSTER")   ~ to,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(lc)) %>%
  distinct(lc) %>%
  pull(lc)

cat("Found", length(lc_gv_connected), "LC contigs connected to GVs\n\n")


# Build full graph --------------------------------------------------------
cat("=== Building graph and calculating centrality ===\n")

g <- graph_from_data_frame(big_connection_df_filtered, directed = FALSE)

cat("Graph has", vcount(g), "nodes and", ecount(g), "edges\n")


# Calculate degree --------------------------------------------------------
cat("Calculating degree centrality...\n")
deg <- degree(g, mode = "all")


# Calculate betweenness (parallelized) -----------------------------------
cat("Calculating betweenness centrality (parallelized, cutoff = 5)...\n")

# Split vertices into chunks for parallel processing
n_vertices <- vcount(g)
vertex_chunks <- split(1:n_vertices, cut(1:n_vertices, breaks = cores_to_use, labels = FALSE))

# Parallel computation
btw_list <- foreach(chunk = vertex_chunks, 
                    .packages = "igraph",
                    .combine = c) %dopar% {
                      estimate_betweenness(g, v = chunk, directed = FALSE, cutoff = 5)
                    }

btw <- btw_list

cat("Centrality calculations complete\n\n")


# Combine into dataframe --------------------------------------------------
deg_df <- enframe(deg, name = "contig_id", value = "degree")
btw_df <- enframe(btw, name = "contig_id", value = "betweeness")

network_df <- left_join(deg_df, btw_df)
network_df$contig_type <- nodes_df$type[match(network_df$contig_id, nodes_df$contig_id)]  

# Update contig_type to distinguish LC-GV-connected from other LCs
network_df <- network_df %>%
  mutate(contig_type = case_when(
    contig_type == "lc" & contig_id %in% lc_gv_connected ~ "lc_gv_connected",  
    TRUE ~ contig_type
  ))

# Set factor order
network_df$contig_type <- factor(network_df$contig_type,
                                 levels = c("lc", "lc_gv_connected", "gv", "plv", "vph"))


# Summary statistics ------------------------------------------------------
cat("=== Centrality summary statistics ===\n")
summary_stats <- network_df %>%
  group_by(contig_type) %>%
  summarise(
    n = n(),
    degree_median = median(degree),
    degree_mean = mean(degree),
    degree_max = max(degree),
    betweenness_median = median(betweeness),
    betweenness_mean = mean(betweeness),
    betweenness_max = max(betweeness)
  )
print(summary_stats)
cat("\n")


# Plot degree distribution ------------------------------------------------
cat("Creating degree distribution plot...\n")

deg_plot <- ggplot(network_df, aes(x = contig_type, y = degree, fill = contig_type)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("lc", "lc_gv_connected"), c("lc", "gv")),
    map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
  scale_y_log10() +
  labs(
    title = "Degree Distribution",
    subtitle = "All edge types combined", 
    x = "Contig Type",
    y = "Degree (log10)"
  ) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("MC", "MC (NCV-connected)", "NCV", "PLV", "VPH")) +
  scale_fill_manual(values = c(
    vph = "goldenrod1",
    lc = "seagreen",
    plv = "hotpink",
    gv = "steelblue",
    lc_gv_connected = "#248191"
  ))


# Plot betweenness distribution -------------------------------------------
cat("Creating betweenness distribution plot...\n")

bet_plot <- ggplot(network_df, aes(x = contig_type, y = betweeness, fill = contig_type)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("lc", "lc_gv_connected"), c("lc", "gv")),
    map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
  scale_y_log10() +
  labs(
    title = "Betweenness Distribution",
    subtitle = "All edge types combined",
    x = "Contig Type",
    y = "Betweeness (log10)"
  ) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("MC", "MC (NCV-connected)", "NCV", "PLV", "VPH")) +
  scale_fill_manual(values = c(
    vph = "goldenrod1",
    lc = "seagreen",
    plv = "hotpink",
    gv = "steelblue",
    lc_gv_connected = "#248191"
  ))


# Save combined plot ------------------------------------------------------
cat("Saving combined plot...\n")

p <- deg_plot + bet_plot

ggsave(plot = p, file = "final/centrality_plots/centrality_all_layers_lc_vs_lc-gv-connected.png", 
       width = 7, height = 4)
ggsave(plot = p, file = "final/centrality_plots/centrality_all_layers_lc_vs_lc-gv-connected.pdf", 
       width = 7, height = 4)
ggsave(plot = p, file = "final/centrality_plots/centrality_all_layers_lc_vs_lc-gv-connected.svg", 
       width = 7, height = 4)


# Per-layer centrality analysis (optional) --------------------------------
cat("\n=== Running per-layer centrality analysis ===\n")

for(layer in unique(big_connection_df_filtered$type)){
  cat("Processing layer:", layer, "\n")
  
  g_layer <- graph_from_data_frame(big_connection_df_filtered %>% filter(type == layer), 
                                   directed = FALSE)
  
  deg_layer <- degree(g_layer, mode = "all")
  btw_layer <- betweenness(g_layer, directed = FALSE, normalized = TRUE)
  
  deg_df_layer <- enframe(deg_layer, name = "contig_id", value = "degree")
  btw_df_layer <- enframe(btw_layer, name = "contig_id", value = "betweeness")
  
  network_df_layer <- left_join(deg_df_layer, btw_df_layer)
  network_df_layer$contig_type <- nodes_df$type[match(network_df_layer$contig_id, nodes_df$contig_id)]
  # Set factor order
  network_df_layer$contig_type <- factor(network_df_layer$contig_type,
                                         levels = c("lc", "gv", "plv", "vph"))
  
  # Degree plot
  deg_plot_layer <- ggplot(network_df_layer, aes(x = contig_type, y = degree, fill = contig_type)) +
    geom_boxplot() +
    geom_signif(
      comparisons = list(c("lc", "vph"), c("lc", "plv"), c("lc", "gv")),
      map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
    scale_y_log10() +
    labs(
      title = "Degree Distribution",
      subtitle = paste0("Layer: ", layer), 
      x = "Contig Type",
      y = "Degree (log10)"
    ) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("MC", "NCV", "PLV", "VPH")) +
    scale_fill_manual(values = c(
      vph = "goldenrod1",
      lc = "seagreen",
      plv = "hotpink",
      gv = "steelblue"
    ))
  
  # Betweenness plot
  bet_plot_layer <- ggplot(network_df_layer, aes(x = contig_type, y = betweeness, fill = contig_type)) +
    geom_boxplot() +
    geom_signif(
      comparisons = list(c("lc", "vph"), c("lc", "plv"), c("lc", "gv")),
      map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
    scale_y_log10() +
    labs(
      title = "Betweenness Distribution",
      subtitle = paste0("Layer: ", layer),
      x = "Contig Type",
      y = "Betweeness (log10)"
    ) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("MC", "NCV", "PLV", "VPH")) +
    scale_fill_manual(values = c(
      vph = "goldenrod1",
      lc = "seagreen",
      plv = "hotpink",
      gv = "steelblue"
    ))
  
  p_layer <- deg_plot_layer + bet_plot_layer
  
  ggsave(plot = p_layer, file = paste0("final/centrality_plots/centrality_", layer, "_plot.png"), 
         width = 7, height = 4)
  ggsave(plot = p_layer, file = paste0("final/centrality_plots/centrality_", layer, "_plot.pdf"), 
         width = 7, height = 4)
  ggsave(plot = p_layer, file = paste0("final/centrality_plots/centrality_", layer, "_plot.svg"), 
         width = 7, height = 4)
}


cat("\n=== Centrality analysis complete ===\n")
cat("Plots saved to: final/centrality_plots/\n")
