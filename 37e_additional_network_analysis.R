# Author: dlu @ veelab
# Version: 2025-06-18
# Script: 37e_additional_network_analyses.R
# Purpose: Additional network analyses including edge type composition, 
#          node connectivity patterns, hubs, assortativity, and network diagnostics

# Packages ----------------------------------------------------------------
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(igraph)
library(cowplot)
library(patchwork)

# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Create output directory
dir.create("testing/additional_analyses", recursive = TRUE, showWarnings = FALSE)


# Load prepared data ------------------------------------------------------
big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")

cat("=== Additional Network Analyses ===\n\n")


# 1. Edge type composition analysis ---------------------------------------
cat("1. EDGE TYPE COMPOSITION\n")
cat("========================\n")

edge_type_summary <- big_connection_df_filtered %>%
  group_by(type) %>%
  summarise(
    n_edges = n(),
    pct = round(n() / nrow(big_connection_df_filtered) * 100, 2)
  ) %>%
  arrange(desc(n_edges))

print(edge_type_summary)

# Plot edge type distribution
edge_type_plot <- ggplot(edge_type_summary, aes(x = reorder(type, -n_edges), y = n_edges, fill = type)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.5, size = 3) +
  scale_y_log10() +
  labs(
    title = "Edge Type Distribution",
    subtitle = paste0("Total edges: ", nrow(big_connection_df_filtered)),
    x = "Edge Type",
    y = "Number of Edges (log10)"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(edge_type_plot, filename = "testing/additional_analyses/edge_type_distribution.png", 
       width = 8, height = 5)
ggsave(edge_type_plot, filename = "testing/additional_analyses/edge_type_distribution.pdf", 
       width = 8, height = 5)


# 2. Node type interaction patterns ---------------------------------------
cat("\n2. NODE TYPE INTERACTION PATTERNS\n")
cat("==================================\n")

interaction_matrix <- big_connection_df_filtered %>%
  mutate(
    interaction = paste0(
      pmin(from_type, to_type), 
      " <-> ", 
      pmax(from_type, to_type)
    )
  ) %>%
  group_by(interaction, type) %>%
  summarise(n = n(), .groups = "drop")

# Create heatmap
interaction_wide <- interaction_matrix %>%
  tidyr::pivot_wider(names_from = type, values_from = n, values_fill = 0)

print(interaction_wide)

# Overall interaction counts (summed across all edge types)
interaction_totals <- big_connection_df_filtered %>%
  mutate(
    interaction = paste0(
      pmin(from_type, to_type), 
      " <-> ", 
      pmax(from_type, to_type)
    )
  ) %>%
  group_by(interaction) %>%
  summarise(total_edges = n()) %>%
  arrange(desc(total_edges))

print(interaction_totals)

interaction_plot <- ggplot(interaction_totals, aes(x = reorder(interaction, -total_edges), 
                                                   y = total_edges)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = total_edges), vjust = -0.5, size = 3) +
  scale_y_log10() +
  labs(
    title = "Node Type Interaction Patterns",
    subtitle = "Which node types connect to which?",
    x = "Interaction Type",
    y = "Number of Edges (log10)"
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(interaction_plot, filename = "testing/additional_analyses/interaction_patterns.png", 
       width = 10, height = 5)
ggsave(interaction_plot, filename = "testing/additional_analyses/interaction_patterns.pdf", 
       width = 10, height = 5)


# 3. Network connectivity diagnostics -------------------------------------
cat("\n3. NETWORK CONNECTIVITY DIAGNOSTICS\n")
cat("====================================\n")

g <- graph_from_data_frame(big_connection_df_filtered, directed = FALSE)

# Basic connectivity metrics
cat("Basic connectivity:\n")
cat("  Connected:", is.connected(g), "\n")
cat("  Number of components:", count_components(g), "\n")

components <- components(g)
component_sizes <- table(components$membership)

cat("  Largest component size:", max(component_sizes), "nodes\n")
cat("  Smallest component size:", min(component_sizes), "node(s)\n")
cat("  Components with >10 nodes:", sum(component_sizes > 10), "\n")
cat("  Isolated nodes (singletons):", sum(component_sizes == 1), "\n")

# Plot component size distribution
component_df <- data.frame(
  component_id = as.numeric(names(component_sizes)),
  size = as.numeric(component_sizes)
) %>%
  arrange(desc(size))

component_plot <- ggplot(component_df, aes(x = 1:nrow(component_df), y = size)) +
  geom_col(fill = "steelblue") +
  scale_y_log10() +
  labs(
    title = "Connected Component Sizes",
    subtitle = paste0(count_components(g), " components total"),
    x = "Component Rank",
    y = "Component Size (log10)"
  ) +
  theme_cowplot()

ggsave(component_plot, filename = "testing/additional_analyses/component_sizes.png", 
       width = 8, height = 5)
ggsave(component_plot, filename = "testing/additional_analyses/component_sizes.pdf", 
       width = 8, height = 5)


# 4. Node degree distribution by type -------------------------------------
cat("\n4. NODE DEGREE DISTRIBUTION BY TYPE\n")
cat("====================================\n")

deg <- degree(g)
deg_df <- data.frame(
  contig_id = names(deg),
  degree = deg
)
deg_df$type <- contig_df$type[match(deg_df$contig_id, contig_df$contig_id)]

deg_summary <- deg_df %>%
  group_by(type) %>%
  summarise(
    n_nodes = n(),
    min_deg = min(degree),
    median_deg = median(degree),
    mean_deg = round(mean(degree), 2),
    max_deg = max(degree),
    n_isolated = sum(degree == 0)
  )

print(deg_summary)

# Violin plot of degree distributions
deg_violin <- ggplot(deg_df %>% filter(degree > 0), 
                     aes(x = type, y = degree, fill = type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  scale_y_log10() +
  scale_fill_manual(values = c(
    vph = "goldenrod1",
    lc = "seagreen",
    plv = "hotpink",
    gv = "steelblue"
  )) +
  labs(
    title = "Degree Distribution by Node Type",
    subtitle = "Violin plots showing full distribution",
    x = "Node Type",
    y = "Degree (log10)"
  ) +
  scale_x_discrete(labels = c("MC", "NCV", "PLV", "VPH")) +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave(deg_violin, filename = "testing/additional_analyses/degree_violin_by_type.png", 
       width = 7, height = 5)
ggsave(deg_violin, filename = "testing/additional_analyses/degree_violin_by_type.pdf", 
       width = 7, height = 5)


# 5. Hub node identification ----------------------------------------------
cat("\n5. HUB NODE IDENTIFICATION\n")
cat("===========================\n")

# Define hubs as nodes with degree > 95th percentile
hub_threshold <- quantile(deg, 0.95)
cat("Hub threshold (95th percentile):", hub_threshold, "\n")

hubs <- deg_df %>%
  filter(degree >= hub_threshold) %>%
  arrange(desc(degree))

cat("Number of hub nodes:", nrow(hubs), "\n")
cat("Hub node types:\n")
print(table(hubs$type))

cat("\nTop 20 hub nodes:\n")
print(head(hubs, 20))

# Save hub nodes to file
fwrite(hubs, "testing/additional_analyses/hub_nodes.csv")
cat("\nHub nodes saved to: testing/additional_analyses/hub_nodes.csv\n")


# 6. Assortativity analysis -----------------------------------------------
cat("\n6. ASSORTATIVITY ANALYSIS\n")
cat("==========================\n")

# Assortativity by node type
V(g)$type_numeric <- as.numeric(factor(contig_df$type[match(V(g)$name, contig_df$contig_id)]))
assortativity_type <- assortativity_nominal(g, V(g)$type_numeric, directed = FALSE)

cat("Assortativity by node type:", round(assortativity_type, 4), "\n")
if(assortativity_type > 0){
  cat("  Interpretation: Nodes tend to connect to similar types (assortative)\n")
} else {
  cat("  Interpretation: Nodes tend to connect to different types (disassortative)\n")
}

# Assortativity by degree
assortativity_deg <- assortativity_degree(g, directed = FALSE)
cat("\nAssortativity by degree:", round(assortativity_deg, 4), "\n")
if(assortativity_deg > 0){
  cat("  Interpretation: High-degree nodes connect to high-degree nodes (assortative)\n")
} else {
  cat("  Interpretation: High-degree nodes connect to low-degree nodes (disassortative)\n")
}


# 7. Network density by edge type -----------------------------------------
cat("\n7. NETWORK DENSITY BY EDGE TYPE\n")
cat("================================\n")

density_by_type <- data.frame(
  edge_type = character(),
  n_edges = numeric(),
  density = numeric(),
  stringsAsFactors = FALSE
)

for(etype in unique(big_connection_df_filtered$type)){
  g_sub <- graph_from_data_frame(
    big_connection_df_filtered %>% filter(type == etype), 
    directed = FALSE
  )
  
  density_by_type <- rbind(density_by_type, data.frame(
    edge_type = etype,
    n_edges = ecount(g_sub),
    density = edge_density(g_sub)
  ))
}

density_by_type <- density_by_type %>% arrange(desc(density))
print(density_by_type)


# Summary report ----------------------------------------------------------
cat("\n=== Additional analyses complete ===\n")
cat("Results saved to: testing/additional_analyses/\n\n")

cat("KEY FINDINGS:\n")
cat("-------------\n")
cat("Total network has", vcount(g), "nodes and", ecount(g), "edges\n")
cat("Network density:", round(edge_density(g), 6), "\n")
cat("Largest component contains", max(component_sizes), "nodes\n")
cat("Node type assortativity:", round(assortativity_type, 4), 
    ifelse(assortativity_type > 0, "(assortative)", "(disassortative)"), "\n")
cat("Degree assortativity:", round(assortativity_deg, 4),
    ifelse(assortativity_deg > 0, "(assortative)", "(disassortative)"), "\n")
cat("\nMost common interactions:\n")
print(head(interaction_totals, 5))

# Clean up
rm(g)
