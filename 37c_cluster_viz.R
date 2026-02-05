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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Create output directories
dir.create("final/plv_vph_subclusters", recursive = TRUE, showWarnings = FALSE)
dir.create("final/gv_subclusters", recursive = TRUE, showWarnings = FALSE)


# Load prepared data ------------------------------------------------------
big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")
GV_info <- fread("intermediate/network/network_analysis/GV_info.csv")
vph_plv_combined_info <- fread("intermediate/network/network_analysis/vph_plv_combined_info.csv")
crispr_df <- fread("intermediate/network/network_analysis/crispr_df.csv")

cat("Loaded data for subgraph visualization\n\n")


# User configuration ------------------------------------------------------
# Set which contig to visualize (single mode) or batch mode
SINGLE_MODE <- FALSE  # Set to TRUE to visualize just one contig
CONTIG_OF_INTEREST <- "Bjer_2_3"  # Only used if SINGLE_MODE = TRUE
SAVE_PLOT <- TRUE

# For batch mode: which sequences to process
plvs <- str_remove(list.files("intermediate/contigs/plv"), "\\.fna")
vphs <- str_remove(list.files("intermediate/contigs/vph"), "\\.fna")
gvs <- GV_info$shortname

# Choose which set to process (change as needed)
list_of_sequences_to_print <- plvs  # Change to: vphs, gvs, or c(plvs, vphs, gvs) for all


# Build graph with Louvain clustering -------------------------------------
cat("Building graph and performing Louvain clustering (resolution = 3)...\n")

graph <- as_tbl_graph(big_connection_df_filtered, directed = FALSE)
graph <- graph %>% 
  mutate(louvain_group = group_louvain(resolution = 3))

cat("Graph has", vcount(graph), "nodes and", ecount(graph), "edges\n")
cat("Found", max(V(graph)$louvain_group), "Louvain communities\n\n")


# Process sequences -------------------------------------------------------
if(SINGLE_MODE){
  list_of_sequences_to_print <- CONTIG_OF_INTEREST
  cat("=== SINGLE MODE: Processing", CONTIG_OF_INTEREST, "===\n\n")
} else {
  cat("=== BATCH MODE: Processing", length(list_of_sequences_to_print), "sequences ===\n\n")
}

n_processed <- 0
n_skipped_no_group <- 0
n_skipped_too_large <- 0

for(contig in list_of_sequences_to_print){
  
  # Get the Louvain group of the current contig
  target_group <- graph %>% 
    as_tibble() %>% 
    filter(name == contig) %>%
    pull(louvain_group)
  
  # If this contig is not in any subcluster, skip it
  if(identical(target_group, integer(0))){
    cat("  SKIP (not in graph):", contig, "\n")
    n_skipped_no_group <- n_skipped_no_group + 1
    next
  }
  
  # Filter the graph to contain only this Louvain group
  subgraph <- graph %>%
    filter(louvain_group == target_group)
  
  # Construct df that contains info on all contigs in this subgraph
  small_node_df <- data.table(id = names(V(subgraph)),
                              contig_type = "")
  
  # If this subgraph is too large (>= 500 nodes), skip
  if(nrow(small_node_df) >= 500){
    cat("  SKIP (too large):", contig, " - ", nrow(small_node_df), "nodes\n")
    n_skipped_too_large <- n_skipped_too_large + 1
    next
  }
  
  # Add contig type
  small_node_df$contig_type <- contig_df$type[match(small_node_df$id, contig_df$contig_id)]
  small_node_df$public_ID <- small_node_df$id
  small_node_df$public_ID[small_node_df$contig_type == "gv"] <- 
    GV_info$public_ID[match(small_node_df$id[small_node_df$contig_type == "gv"], GV_info$shortname)]
  
  # Add CRISPR info
  small_node_df$cas_genes <- ""
  for(i in 1:nrow(small_node_df)){
    # Virophages get an NA
    if(small_node_df$contig_type[i] == "vph"){
      small_node_df$cas_genes[i] <- "Not applicable"
      next
    }
    
    # GVs and LCs get a concatenated unique list of cas genes found
    tmp_df <- crispr_df %>% 
      filter(contig_id == small_node_df$id[i])
    
    if(nrow(tmp_df) >= 1){
      cas_genes <- tmp_df %>% 
        select(annotation) %>% 
        unlist() %>% 
        unique() %>% 
        paste(collapse = ", ")
    } else {
      cas_genes <- "None detected."
    }
    small_node_df$cas_genes[i] <- cas_genes
  }
  
  # Subset edges to only those in this community
  small_edge_df <- big_connection_df_filtered %>%
    filter(from %in% small_node_df$id & to %in% small_node_df$id)
  
  # Build subgraph
  g <- graph_from_data_frame(small_edge_df, directed = FALSE, small_node_df)
  V(g)$node_type <- contig_df$type[match(V(g)$name, contig_df$contig_id)]
  V(g)$tax_info <- contig_df$tax_info[match(V(g)$name, contig_df$contig_id)]
  E(g)$edge_type <- as.factor(E(g)$type)
  
  g <- as_tbl_graph(g)
  
  # Highlight the focal contig
  g <- g %>%
    activate(nodes) %>%
    mutate(is_focal = (name == contig))
  
  # Create layout
  layout <- create_layout(g, layout = "fr")
  
  # Plot
  ggiraph_plot <- ggraph(layout) +
    geom_edge_link(aes(color = type), show.legend = FALSE) +
    geom_point_interactive(shape = 21,
                           aes(x = x, 
                               y = y,
                               size = if_else(is_focal, 4, 4),
                               alpha = if_else(is_focal, 1, 0.8),   
                               fill = node_type,
                               data_id = public_ID,
                               tooltip = paste("Node ID:", public_ID,
                                               "\nType:", node_type,
                                               "\nTaxonomy:", tax_info,
                                               "\nCas genes:", cas_genes),
                               color = if_else(is_focal, "black", "white"),
                               stroke = if_else(is_focal, 1, 0)
                           )) +
    scale_fill_manual(values = c(
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
    name = "Node Type") +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_color_identity() +
    facet_edges(~edge_type, ncol = 3, labeller = labeller(
      edge_type = c(
        "gene_sharing" = "Gene Sharing",
        "crispr" = "CRISPR",
        "non_crispr" = "MIMIVIRE",
        "integration_boundary" = "Integration (boundary)",
        "integration_middle" = "Integration (middle)",
        "occurance_positive" = "Co-occurrence (positive)", 
        "occurance_negative" = "Co-occurrence (negative)"
      )
    )) +
    theme_classic() +
    theme(legend.position = "bottom") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
  
  interactive_plot <- girafe(ggobj = ggiraph_plot)
  
  # Save plots
  if(SAVE_PLOT){
    type_of_contig <- str_remove(contig, "^.*\\_")
    
    if(type_of_contig == "plv" || type_of_contig == "vph"){
      public_ID <- vph_plv_combined_info$public_ID[vph_plv_combined_info$contig_ID == contig]
      FILE_BASE <- paste0("final/plv_vph_subclusters/", public_ID) 
    } else {
      public_ID <- GV_info$public_ID[GV_info$shortname == contig]
      FILE_BASE <- paste0("final/gv_subclusters/", public_ID)
    }
    
    # Adjust height based on number of facets
    FACETS <- length(unique(small_edge_df$type))
    HEIGHT_FACTOR <- ceiling(FACETS / 3)
    
    ggsave(ggiraph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".svg"), 
           width = 4.5, height = 1.8 * HEIGHT_FACTOR)
    ggsave(ggiraph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".pdf"), 
           width = 4.5, height = 1.8 * HEIGHT_FACTOR)
    ggsave(ggiraph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".png"), 
           width = 4.5, height = 1.8 * HEIGHT_FACTOR)
    htmltools::save_html(interactive_plot, file = paste0(FILE_BASE, "_interactive.html"))
    
    cat("  PROCESSED:", contig, "(", nrow(small_node_df), "nodes,", 
        nrow(small_edge_df), "edges )\n")
    n_processed <- n_processed + 1
  }
}

# Summary
cat("\n=== Subgraph visualization complete ===\n")
cat("Processed:", n_processed, "contigs\n")
cat("Skipped (not in graph):", n_skipped_no_group, "\n")
cat("Skipped (too large):", n_skipped_too_large, "\n")
cat("Plots saved to: final/plv_vph_subclusters/ and final/gv_subclusters/\n")
