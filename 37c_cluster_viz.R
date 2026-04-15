# Author: dlu @ veelab
# Version: 2025-06-18
# Script: 37c_visualize_node_subgraphs.R
# Purpose: Visualize network clusters surrounding specific nodes (PLVs, VPHs, GVs)

# Packages ----------------------------------------------------------------
library(dplyr)
library(data.table)
library(stringr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggiraph)
library(htmltools)

# Setup -------------------------------------------------------------------
setwd("/run/user/1000/gvfs/sftp:host=login01.lisc.univie.ac.at,user=luecking/lisc/home/user/luecking/luecking_scratch/projects/wwtp_linking/")
rm(list = ls())
# dir.create("final/plv_vph_subclusters", recursive = TRUE, showWarnings = FALSE)
# dir.create("final/gv_subclusters", recursive = TRUE, showWarnings = FALSE)


# Load prepared data ------------------------------------------------------
big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")
GV_info <- fread("intermediate/network/network_analysis/GV_info.csv")
vph_plv_combined_info <- fread("intermediate/network/network_analysis/vph_plv_combined_info.csv")
crispr_df <- fread("intermediate/network/network_analysis/crispr_df.csv")
lc_tax_info_short <- fread("intermediate/network/network_analysis/lc_tax_info_short.csv")
membership_df <- fread("intermediate/mash/lc_contigs_only_cluster_membership_005.csv")

membership_df$contig_id <- str_remove(membership_df$contig_id, "\\spolypolish")
membership_df$cluster_id <- str_remove(membership_df$cluster_id, "\\spolypolish")

cat("Loaded data for subgraph visualization\n\n")


# Create cluster-level taxonomy -------------------------------------------
cat("Creating cluster-level taxonomy...\n")

cluster_taxonomy <- membership_df %>%
  filter(grepl("LC_CLUSTER", cluster_id)) %>%
  left_join(lc_tax_info_short %>% select(contig_id, majority_organism), by = "contig_id") %>%
  filter(!is.na(majority_organism) & majority_organism != "") %>%
  group_by(cluster_id) %>%
  summarize(
    cluster_size = n(),
    tax_list = list(majority_organism),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    tax_info = {
      tax_vec <- unlist(tax_list)
      counts <- sort(table(tax_vec), decreasing = TRUE)
      total <- sum(counts)
      
      if(total <= 3) {
        paste(paste0(names(counts), " (", counts, ")"), collapse = "; ")
      } else {
        counts_filtered <- counts[counts >= 3]
        if(length(counts_filtered) > 0) {
          paste(paste0(names(counts_filtered), " (", counts_filtered, ")"), collapse = "; ")
        } else {
          "Mixed <3 each"
        }
      }
    }
  ) %>%
  ungroup() %>%
  select(cluster_id, cluster_size, tax_info)

cat("Created taxonomy for", nrow(cluster_taxonomy), "LC clusters\n\n")

# Create nodes_df for graph construction ----------------------------------
nodes_df <- contig_df %>%
  left_join(membership_df %>% select(contig_id, cluster_id, cluster_size), by = "contig_id") %>%
  mutate(
    node_id = ifelse(!is.na(cluster_id), cluster_id, contig_id),
    original_contig_id = contig_id
  ) %>%
  as_tibble() %>%
  distinct(node_id, .keep_all = TRUE) %>%
  select(-contig_id) %>%
  select(contig_id = node_id, everything()) %>%  # Combined rename + select
  mutate(
    type = case_when(
      grepl("LC_CLUSTER", contig_id) ~ "lc",
      str_ends(contig_id, "vph") ~ "vph",
      str_ends(original_contig_id, "lc") ~ "lc",
      str_ends(contig_id, "plv") ~ "plv",
      TRUE ~ "gv"
    )
  )

# Update taxonomy for LC clusters
nodes_df <- nodes_df %>%
  left_join(cluster_taxonomy %>% select(cluster_id, tax_info_cluster = tax_info),
            by = c("contig_id" = "cluster_id")) %>%
  mutate(tax_info = ifelse(!is.na(tax_info_cluster), tax_info_cluster, tax_info)) %>%
  select(-tax_info_cluster)


# Build graph with Louvain clustering -------------------------------------
cat("Building graph and performing Louvain clustering (resolution = 10)...\n")

graph <- graph_from_data_frame(
  d = big_connection_df_filtered %>% select(from, to, type, weight),
  directed = FALSE,
  vertices = nodes_df
)

graph <- as_tbl_graph(graph)
set.seed(1312) # I need to set a seed, otherwise this is non-deterministic
graph <- graph %>% 
  mutate(louvain_group = group_louvain(weights = weight, resolution = 10))

cat("Graph has", vcount(graph), "nodes and", ecount(graph), "edges\n")
cat("Found", max(V(graph)$louvain_group), "Louvain communities\n\n")


# User configuration ------------------------------------------------------
SINGLE_MODE <- FALSE
CONTIG_OF_INTEREST <- "Bjer_2_3"
SAVE_PLOT <- TRUE

plvs <- str_remove(list.files("intermediate/contigs/plv"), "\\.fna")
vphs <- str_remove(list.files("intermediate/contigs/vph"), "\\.fna")
plv_vph <- c(plvs, vphs)

gvs <- GV_info$shortname

list_of_sequences_to_print <- plv_vph


# Process sequences -------------------------------------------------------
if(SINGLE_MODE){
  list_of_sequences_to_print <- CONTIG_OF_INTEREST
  contig <- CONTIG_OF_INTEREST
  cat("=== SINGLE MODE: Processing", CONTIG_OF_INTEREST, "===\n\n")
} else {
  cat("=== BATCH MODE: Processing", length(list_of_sequences_to_print), "sequences ===\n\n")
}

n_processed <- 0
n_skipped_no_group <- 0
n_skipped_too_large <- 0

for(contig in list_of_sequences_to_print){
  
  target_group <- graph %>% 
    as_tibble() %>% 
    filter(name == contig) %>%
    pull(louvain_group)
  
  if(identical(target_group, integer(0))){
    cat("  SKIP (not in graph):", contig, "\n")
    n_skipped_no_group <- n_skipped_no_group + 1
    next
  }
  
  subgraph <- graph %>%
    filter(louvain_group == target_group)
  
  small_node_df <- data.table(id = names(V(subgraph)), contig_type = "")
  
  if(nrow(small_node_df) >= 500){
    cat("  SKIP (too large):", contig, " - ", nrow(small_node_df), "nodes\n")
    n_skipped_too_large <- n_skipped_too_large + 1
    next
  }
  
  # Get type from nodes_df (which has clusters)
  small_node_df$contig_type <- nodes_df$type[match(small_node_df$id, nodes_df$contig_id)]
  
  # Get public IDs with NA handling
  small_node_df$public_ID <- small_node_df$id
  
  # Add CRISPR info
  small_node_df$cas_genes <- ""
  for(i in 1:nrow(small_node_df)){
    if(small_node_df$contig_type[i] == "vph"){
      small_node_df$cas_genes[i] <- "Not applicable"
      next
    }
    
    # For LC clusters, check all members
    if(grepl("LC_CLUSTER", small_node_df$id[i])){
      member_ids <- membership_df %>% 
        filter(cluster_id == small_node_df$id[i]) %>% 
        pull(contig_id)
      tmp_df <- crispr_df %>% filter(contig_id %in% member_ids)
    } else {
      tmp_df <- crispr_df %>% filter(contig_id == small_node_df$id[i])
    }
    
    if(nrow(tmp_df) >= 1){
      cas_genes <- tmp_df %>% 
        pull(annotation) %>% 
        unique() %>% 
        paste(collapse = ", ")
    } else {
      cas_genes <- "None detected"
    }
    small_node_df$cas_genes[i] <- cas_genes
  }
  
  small_edge_df <- big_connection_df_filtered %>%
    filter(from %in% small_node_df$id & to %in% small_node_df$id)
  
  # Skip if no edges
  if(nrow(small_edge_df) == 0){
    cat("  SKIP (no edges):", contig, "\n")
    n_skipped_no_group <- n_skipped_no_group + 1
    next
  }
  
  # Add is_directed flag to edges
  directed_types <- c("crispr", "integration_boundary", "integration_middle")
  small_edge_df$is_directed <- small_edge_df$type %in% directed_types
  
  # Build subgraph (undirected)
  g <- graph_from_data_frame(small_edge_df, 
                             directed = TRUE, 
                             small_node_df)
  
  V(g)$node_type <- nodes_df$type[match(V(g)$name, nodes_df$contig_id)]
  V(g)$tax_info <- nodes_df$tax_info[match(V(g)$name, nodes_df$contig_id)]
  E(g)$edge_type <- as.factor(E(g)$type)
  
  g <- as_tbl_graph(g)
  g <- g %>%
    activate(nodes) %>%
    mutate(is_focal = (name == contig)) %>%
    activate(edges) %>%  # Switch to edges!
    mutate(arrow_size = if_else(type %in% directed_types, 2, 0))
  
  layout <- create_layout(g, layout = "gem")
  
  graph_plot <- ggraph(layout) +
    # Layer 1: Undirected edges (no filter, no arrow)
    geom_edge_link0(
      aes(color = type, filter = arrow_size == 0),
      alpha = 0.9,
      show.legend = FALSE
    ) +
    # Layer 2: Directed edges (with arrow)
    geom_edge_link0(
      aes(color = type, filter = arrow_size > 0),
      arrow = arrow(length = unit(2, 'mm'), type = "closed"),
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_point_interactive(
      shape = 21,
      aes(
        x = x, 
        y = y,
        size = if_else(is_focal, 4, 3),
        alpha = if_else(is_focal, 1, 0.8),   
        fill = node_type,
        data_id = public_ID,
        tooltip = paste("Node ID:", public_ID,
                        "\nType:", node_type,
                        "\nTaxonomy:", tax_info,
                        "\nCas genes:", cas_genes),
        color = if_else(is_focal, "black", "black"),
        stroke = if_else(is_focal, 1, 0.3)
      )
    ) +
    scale_fill_manual(
      values = c(
        vph = "goldenrod1",
        lc = "seagreen",
        plv = "hotpink",
        gv = "steelblue"
      ),
      labels = c(
        vph = "VPH",
        lc = "MC",
        plv = "PLV",
        gv = "NCV"
      ),
      name = "Node Type"
    ) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_color_identity() +
    facet_edges(
      ~edge_type,
      ncol = 3,
      labeller = labeller(
        edge_type = c(
          "gene_sharing" = "Gene Sharing",
          "crispr" = "CRISPR",
          # "non_crispr" = "Repeats", # removed!
          "integration_boundary" = "Integration (boundary)",
          "integration_middle" = "Integration (middle)",
          "occurance_positive" = "Co-occurrence (positive)", 
          "occurance_negative" = "Co-occurrence (negative)"
        )
      )
    ) +
    theme_bw() +
    theme(legend.position = "bottom") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  interactive_plot <- girafe(ggobj = graph_plot)
  interactive_plot
  
  if(SAVE_PLOT){
    type_of_contig <- str_remove(contig, "^.*\\_")
    
    if(type_of_contig %in% c("plv", "vph")){
      public_ID <- vph_plv_combined_info$public_ID[vph_plv_combined_info$contig_ID == contig]
      FILE_BASE <- paste0("final/plv_vph_subcluster/", public_ID) 
    } else {
      public_ID <- GV_info$public_ID[GV_info$shortname == contig]
      FILE_BASE <- paste0("final/gv_subcluster/", public_ID)
    }
    
    FACETS <- length(unique(small_edge_df$type))
    HEIGHT_FACTOR <- ceiling(FACETS / 3)
    
    ggsave(graph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".svg"), 
           width = 6, height = 2 * HEIGHT_FACTOR)
    ggsave(graph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".pdf"), 
           width = 6, height = 2 * HEIGHT_FACTOR)
    ggsave(graph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".png"), 
           width = 6, height = 2 * HEIGHT_FACTOR)
    htmltools::save_html(interactive_plot, file = paste0(FILE_BASE, "_interactive.html"))
    
    cat("  PROCESSED:", contig, "(", nrow(small_node_df), "nodes,", 
        nrow(small_edge_df), "edges )\n")
    n_processed <- n_processed + 1
  }
}

cat("\n=== Subgraph visualization complete ===\n")
cat("Processed:", n_processed, "contigs\n")
cat("Skipped (not in graph):", n_skipped_no_group, "\n")
cat("Skipped (too large):", n_skipped_too_large, "\n")
cat("Plots saved to: final/plv_vph_subclusters/ and final/gv_subclusters/\n")
