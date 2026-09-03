# Author: dlu @ veelab
# Version: 2025-06-18
# Script: 37b_community_and_scalefree_analysis.R
# Purpose: Analyze community structure (Louvain) and test scale-free properties

# Packages ----------------------------------------------------------------
library(dplyr)
library(data.table)
library(ggplot2)
library(igraph)
library(tidygraph)
library(poweRlaw)

# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")
rm(list = ls())
# Create output directories
dir.create("intermediate/network/network_analysis/louvain_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("intermediate/network/network_analysis/scale_free_analysis", recursive = TRUE, showWarnings = FALSE)


# Load prepared data ------------------------------------------------------
big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered_REVIEW.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")
membership_df <- fread("intermediate/mash/lc_contigs_only_cluster_membership_005.csv")


# Create cluster-level node dataframe ------------------------------------
cat("Creating cluster-level node dataframe...\n")

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


# =========================================================================
# PART 1: COMMUNITY STRUCTURE ANALYSIS
# =========================================================================

cat("=== PART 1: COMMUNITY STRUCTURE ANALYSIS ===\n")
cat("Testing weighted Louvain clustering with each edge type removed...\n\n")

res <- list()

for(t in c(unique(big_connection_df_filtered$type), "none")){
  cat("Processing:", t, "\n")
  
  if(t != "none"){
    edges_to_use <- big_connection_df_filtered %>% filter(type != t)
  } else {
    edges_to_use <- big_connection_df_filtered
  }
  
  # Build graph from data frame (includes weights)
  graph <- graph_from_data_frame(
    d = edges_to_use, 
    directed = FALSE, 
    vertices = nodes_df
  )
  
  # Weighted Louvain clustering
  communities <- cluster_louvain(graph, weights = E(graph)$weight, resolution = 10)
  
  # Extract membership
  membership_df <- data.frame(
    contig_id = V(graph)$name,
    louvain_group = membership(communities),
    stringsAsFactors = FALSE
  )
  
  # Create summary table
  a <- as.data.table(table(membership_df$louvain_group))
  setnames(a, c("V1", "N"), c("group", "count"))
  a[, group := as.integer(as.character(group))]
  a[, removed_type := t]
  res[[t]] <- a
}

dt <- rbindlist(res)

cat("\n=== Community size summary ===\n")
print(dt[, .(
  n_communities = .N,
  min_size = as.numeric(min(count)),
  median_size = as.numeric(median(count)),
  max_size = as.numeric(max(count)),
  n_singletons = as.numeric(sum(count == 1)),
  n_10_to_500 = as.numeric(sum(count >= 10 & count <= 500))
), by = removed_type])


# Plot ECDF of community sizes --------------------------------------------
cat("\nCreating ECDF plot...\n")

ecdf_plot <- ggplot(dt, aes(x = count, colour = removed_type)) +
  stat_ecdf() +
  scale_x_log10() +
  theme_minimal() +
  labs(
    x = "Community size (log scale)",
    y = "ECDF",
    title = "Distribution of Louvain community sizes",
    subtitle = "How does edge type removal affect community structure?"
  ) +
  theme(legend.position = "bottom")

ggsave(plot = ecdf_plot, filename = "intermediate/network/network_analysis/louvain_plots/ecdf.png", width = 8, height = 6)
ggsave(plot = ecdf_plot, filename = "intermediate/network/network_analysis/louvain_plots/ecdf.pdf", width = 8, height = 6)


# Plot multipanel community size distributions ----------------------------
cat("Creating multipanel plot...\n")

multi_plot <- ggplot(dt, aes(x = group, y = count)) +
  geom_col() +
  facet_wrap(~ removed_type, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Louvain group",
    y = "Number of contigs",
    title = "Louvain community size distributions",
    subtitle = "Panel = this edge type has been removed"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(plot = multi_plot, filename = "intermediate/network/network_analysis/louvain_plots/multipanel.png", width = 12, height = 8)
ggsave(plot = multi_plot, filename = "intermediate/network/network_analysis/louvain_plots/multipanel.pdf", width = 12, height = 8)


# Test different Louvain resolutions --------------------------------------
cat("\n=== Testing different Louvain resolutions ===\n")

# Build full weighted graph for resolution testing
graph <- graph_from_data_frame(
  d = big_connection_df_filtered, 
  directed = FALSE, 
  vertices = nodes_df
)

for(res_val in seq(0.5, 15, by = 1)) {
  communities <- cluster_louvain(graph, weights = E(graph)$weight, resolution = res_val)
  sizes <- table(membership(communities))
  
  cat("\nResolution:", res_val, "\n")
  cat("  N communities:", length(sizes), "\n")
  cat("  Largest:", max(sizes), "nodes\n")
  cat("  Communities 10-500 nodes:", sum(sizes >= 10 & sizes <= 500), "\n")
  cat("  Singletons:", sum(sizes == 1), "\n")
}



# =========================================================================
# PART 2: SCALE-FREE ANALYSIS
# =========================================================================

cat("\n\n=== PART 2: SCALE-FREE NETWORK ANALYSIS ===\n")
cat("Building graph and calculating degree distribution...\n")

g <- graph_from_data_frame(d = big_connection_df_filtered, 
                           directed = FALSE, 
                           vertices = nodes_df)

d <- degree(g)

cat("Degree statistics:\n")
cat("  Min:", min(d), "\n")
cat("  Median:", median(d), "\n")
cat("  Mean:", round(mean(d), 2), "\n")
cat("  Max:", max(d), "\n\n")


# Create cumulative degree distribution -----------------------------------
cat("Creating cumulative degree distribution plot...\n")

# Create frequency table and calculate cumulative probability
degree_counts <- as.data.frame(table(d))
names(degree_counts) <- c("k", "Freq")
degree_counts$k <- as.numeric(as.character(degree_counts$k))

# Sort and calculate cumulative distribution P(X >= k)
degree_counts <- degree_counts[order(degree_counts$k, decreasing = TRUE), ]
degree_counts$cumulative_prob <- cumsum(degree_counts$Freq) / sum(degree_counts$Freq)

# Log-log plot to check for scale-free properties
loglog_plot <- ggplot(degree_counts, aes(x = k, y = cumulative_prob)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  scale_x_log10() + 
  scale_y_log10() +
  annotation_logticks() +
  labs(
    title = "Log-Log Degree Distribution",
    subtitle = "Check for linearity to identify scale-free properties",
    x = "Degree (k)",
    y = "Cumulative Probability P(k)"
  ) +
  theme_minimal()

ggsave(plot = loglog_plot, filename = "intermediate/network/network_analysis/scale_free_analysis/loglog_degree_distribution.png", 
       width = 8, height = 6)
ggsave(plot = loglog_plot, filename = "intermediate/network/network_analysis/scale_free_analysis/loglog_degree_distribution.pdf", 
       width = 8, height = 6)


# Test power law distribution with poweRlaw -------------------------------
cat("\n=== Testing power law distribution ===\n")

# Remove zero-degree nodes
d_positive <- d[d > 0]

# Create power law distribution object
m_pl <- displ$new(d_positive)

# Estimate the minimum degree (xmin) where the power law starts
cat("Estimating xmin (minimum degree for power law tail)...\n")
est <- estimate_xmin(m_pl)
m_pl$setXmin(est)

cat("Power law parameters:\n")
cat("  xmin:", m_pl$getXmin(), "\n")
cat("  alpha (scaling exponent):", round(m_pl$getPars(), 3), "\n\n")


# Compare power law to log-normal distribution ---------------------------
cat("Comparing power law to log-normal distribution...\n")

m_ln <- dislnorm$new(d_positive)
m_ln$setXmin(m_pl$getXmin())  # Compare on the same tail
m_ln$setPars(estimate_pars(m_ln))

comp <- compare_distributions(m_pl, m_ln)

cat("\nLikelihood ratio test:\n")
cat("  Test statistic:", round(comp$test_statistic, 4), "\n")
if(comp$test_statistic > 0){
  cat("  Interpretation: Power law is a BETTER fit than log-normal\n")
} else {
  cat("  Interpretation: Log-normal is a BETTER fit than power law\n")
}
cat("  p-value:", round(comp$p_two_sided, 4), "\n")

if(comp$p_two_sided < 0.05){
  cat("  Conclusion: Significant difference between models (p < 0.05)\n")
} else {
  cat("  Conclusion: No significant difference between models (p >= 0.05)\n")
}


# Create comparison plot --------------------------------------------------
cat("\nCreating comparison plot of fitted distributions...\n")

# Generate comparison data
x_vals <- seq(m_pl$getXmin(), max(d_positive), length.out = 100)

# Power law predictions
pl_probs <- dist_cdf(m_pl, x_vals, lower_tail = FALSE)

# Log-normal predictions
ln_probs <- dist_cdf(m_ln, x_vals, lower_tail = FALSE)

# Combine into dataframe
fit_df <- data.frame(
  k = rep(x_vals, 2),
  cumulative_prob = c(pl_probs, ln_probs),
  model = rep(c("Power Law", "Log-Normal"), each = length(x_vals))
)

comparison_plot <- ggplot() +
  geom_point(data = degree_counts %>% filter(k >= m_pl$getXmin()), 
             aes(x = k, y = cumulative_prob), 
             alpha = 0.6, color = "steelblue", size = 2) +
  geom_line(data = fit_df, 
            aes(x = k, y = cumulative_prob, color = model), 
            size = 1) +
  scale_x_log10() + 
  scale_y_log10() +
  annotation_logticks() +
  scale_color_manual(values = c("Power Law" = "red", "Log-Normal" = "darkgreen")) +
  labs(
    title = "Power Law vs Log-Normal Fit",
    subtitle = paste0("Test statistic = ", round(comp$test_statistic, 3), 
                      ", p = ", round(comp$p_two_sided, 4)),
    x = "Degree (k)",
    y = "Cumulative Probability P(k)",
    color = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(plot = comparison_plot, filename = "intermediate/network/network_analysis/scale_free_analysis/power_law_comparison.png", 
       width = 8, height = 6)
ggsave(plot = comparison_plot, filename = "intermediate/network/network_analysis/scale_free_analysis/power_law_comparison.pdf", 
       width = 8, height = 6)


# Final summary -----------------------------------------------------------
cat("\n=== Analysis complete ===\n")
cat("Community structure plots saved to: intermediate/network/network_analysis/scale_free_analysis/\n")
cat("Scale-free analysis saved to: intermediate/network/network_analysis/scale_free_analysis/\n\n")

cat("SUMMARY:\n")
cat("--------\n")
cat("Network has", vcount(g), "nodes and", ecount(g), "edges\n")
cat("Power law tail starts at degree k >=", m_pl$getXmin(), "\n")
cat("Power law exponent (alpha):", round(m_pl$getPars(), 3), "\n")

if(comp$test_statistic > 0 & comp$p_two_sided < 0.05){
  cat("\nNetwork shows STRONG evidence of scale-free properties\n")
} else if(comp$test_statistic > 0){
  cat("\nNetwork shows WEAK evidence of scale-free properties\n")
} else {
  cat("\nNetwork is better described by log-normal than power law\n")
}

# Clean up
rm(g, d, d_positive, m_pl, m_ln, comp, graph)
