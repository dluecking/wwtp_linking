# Author: dlu @ veelab
# Version: 2025-06-18

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(igraph)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# This script collects the information of all layers, which are connecting 
# entities to long_contigs. We need to know, which lcs are connected to 
# something of interest (GV, VPH, PLV), so that we can reduce the occurance 
# layer to these contigs of interest


# load egdgelists of layers that connect to lc ----------------------------

edgelist_crispr <- fread("intermediate/network/crispr.csv") %>% 
  select(from, to)
edgelist_integration_b <- fread("intermediate/network/integration_b.csv")%>% 
  select(from, to)
edgelist_integration_m <- fread("intermediate/network/integration_m.csv")%>% 
  select(from, to)
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing_only_interesting.csv")%>% 
  select(from, to)

big_edgelist <- rbind(edgelist_crispr, edgelist_integration_b, edgelist_integration_m, edgelist_gene_sharing)




# test --------------------------------------------------------------------

# Ensure no duplicate edges and remove self-loops if they exist
g <- graph_from_data_frame(d = unique(big_edgelist), directed = FALSE)

# 3. Assign 'type' attribute to nodes (lc or non-lc)
V(g)$type <- ifelse(grepl("_lc$", V(g)$name), "lc", "non-lc")

# Get list of non-lc contig names
non_lc_nodes <- V(g)$name[V(g)$type == "non-lc"]

# Initialize sets for results
directly_connected_lcs <- character(0)
all_connected_lcs <- character(0) # All LCs connected to any non-lc

# Find all LCs connected to any non-LC node (directly or indirectly)
# Iterate through each non-lc node and find all its neighbors, then their neighbors, etc.
# Keep track of only the 'lc' nodes in these connected components.
for (node_name in non_lc_nodes) {
  # Get the connected component of the current non-lc node
  # `subcomponent` finds all nodes reachable from a given node
  component_nodes <- names(subcomponent(g, node_name, mode = "all"))
  
  # Filter for 'lc' nodes within this component and add to the set
  all_connected_lcs <- union(all_connected_lcs,
                             component_nodes[grepl("_lc$", component_nodes)])
}

# Identify directly connected LCs
# Iterate through edges and check if one end is 'lc' and the other 'non-lc'
edges_list <- as_data_frame(g)
for (i in 1:nrow(edges_list)) {
  node1 <- edges_list$from[i]
  node2 <- edges_list$to[i]
  
  type1 <- V(g)[node1]$type
  type2 <- V(g)[node2]$type
  
  if (type1 == "lc" && type2 == "non-lc") {
    directly_connected_lcs <- union(directly_connected_lcs, node1)
  } else if (type1 == "non-lc" && type2 == "lc") {
    directly_connected_lcs <- union(directly_connected_lcs, node2)
  }
}

# Calculate indirectly connected LCs
# These are LCs that are connected to a non-LC node, but not directly
indirectly_connected_lcs <- setdiff(all_connected_lcs, directly_connected_lcs)

# --- Results ---
message("List of directly connected 'lc' contigs:")
print(directly_connected_lcs)
message("\nNumber of directly connected 'lc' contigs:")
print(length(directly_connected_lcs))

message("\nList of all 'lc' contigs connected (directly or indirectly) to a non-lc contig:")
print(all_connected_lcs)
message("\nTotal number of 'lc' contigs connected (directly or indirectly) to a non-lc contig:")
print(length(all_connected_lcs))

message("\nList of indirectly connected 'lc' contigs:")
print(indirectly_connected_lcs)
message("\nNumber of indirectly connected 'lc' contigs:")
print(length(indirectly_connected_lcs))



# write to file -----------------------------------------------------------

write(all_connected_lcs, "intermediate/network/list_of_connected_lcs.txt")
