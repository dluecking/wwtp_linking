# Author: dlu @ veelab
# Version: 2025-06-18

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(igraph)
library(igraphdata)
library(ggraph)
library(cowplot)
library(tidygraph)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load extra data ---------------------------------------------------------

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

# now prepare a color
contig_df <- contig_df %>% 
  mutate(
    color = case_when(
      str_ends(contig_id, "vph") ~ "goldenrod1",
      str_ends(contig_id, "lc")  ~ "seagreen",
      str_ends(contig_id, "plv") ~ "hotpink",
      TRUE ~ "steelblue"
    )
  )

# load data ---------------------------------------------------------------

edgelist_crispr <- fread("intermediate/network/crispr.csv")
edgelist_integration_b <- fread("intermediate/network/integration_b.csv")
edgelist_integration_m <- fread("intermediate/network/integration_m.csv")
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing_only_interesting.csv")
edgelist_occurance_ill <- fread("intermediate/network/occurance_ill.csv")
edgelist_occurance_ont <- fread("intermediate/network/occurance_ont.csv")
edgelist_non_crispr <- fread("intermediate/network/non_CRISPR.csv")






# explore -----------------------------------------------------------------
# create binary yes/no big list with connection type for each edge
big_connection_df <- rbind(
  edgelist_crispr %>% select(from, to) %>% 
    mutate(type = "crispr"),
  edgelist_gene_sharing %>% select(from, to) %>% 
    mutate(type = "gene_sharing"),
  edgelist_integration_b %>% select(from, to) %>% 
    mutate(type = "integration_boundary"),
  edgelist_integration_m %>% select(from, to) %>% 
    mutate(type = "integration_middle"),
  edgelist_non_crispr %>% select(from, to) %>% 
    mutate(type = "non_crispr"),
  edgelist_occurance_ill %>% filter(correlation == "positive") %>% select(from, to) %>% 
    mutate(type = "occurance_ill_positive"),
  edgelist_occurance_ont %>% filter(correlation == "positive") %>% select(from, to) %>% 
    mutate(type = "occurance_ont_positive"),
  edgelist_occurance_ill %>% filter(correlation == "negative") %>% select(from, to) %>% 
    mutate(type = "occurance_ill_negative"),
  edgelist_occurance_ont %>% filter(correlation == "negative") %>% select(from, to) %>% 
    mutate(type = "occurance_ont_negative")
)




# try to get to facets ----------------------------------------------------

# get a list of ids that are part of the specific subcluster


# this first retains only connections which are of a specific type (e.g. "gene_sharing")
df <- big_connection_df %>% filter(!str_detect(type, ".*occurance.*")) %>% 
  filter(type == "gene_sharing")

graph <- as_tbl_graph(df, directed = F)
graph <- graph %>% 
  mutate(louvain_group = group_louvain())

target_group <- graph %>% 
  as_tibble() %>% 
  filter(name == "Aved_tig00303955-10-54120_vph") %>%
  # filter(name == "Hjor_1") %>%
  pull(louvain_group)
subgraph <- graph %>%
  filter(louvain_group == target_group)

small_node_df <- data.table(id = names(V(subgraph)),
                            contig_type = "")
small_node_df$contig_type <- contig_df$type[match(small_node_df$id, contig_df$contig_id)]

# then subset the big df, so we retain only edges of nodes that are part of the 
# group
small_edge_df <- big_connection_df %>%
  filter(from %in% small_node_df$id & to %in% small_node_df$id)

g <- graph_from_data_frame(small_edge_df, directed = F, small_node_df)
V(g)$node_type <- contig_df$type[match(V(g)$name, contig_df$contig_id)]
E(g)$edge_type <- as.factor(E(g)$type)

g <- as_tbl_graph(g)


layout <- create_layout(g, layout = "fr")
ggraph(layout) +
  geom_edge_link(aes(color = type)) +
  geom_node_point(aes(color = node_type), size = 5, alpha = 0.8) +
  scale_color_manual(values = c(
    vph = "goldenrod1",
    lc = "seagreen",
    plv = "hotpink",
    gv = "steelblue"
  )) +
  facet_edges(~edge_type) 


# another idea: GV-LC partners --------------------------------------------

multi_connection_pairs <- big_connection_df %>%
  group_by(from, to) %>%
  summarise(
    unique_connection_types = n_distinct(type)
  ) %>%
  ungroup() %>% # Ungroup to remove grouping structure
  filter(unique_connection_types > 2)


# now add information about type
multi_connection_pairs$from_type <- contig_df$type[match(multi_connection_pairs$from, contig_df$contig_id)]
multi_connection_pairs$to_type <- contig_df$type[match(multi_connection_pairs$to, contig_df$contig_id)]


multi_connection_pairs$connection_type <- paste0(multi_connection_pairs$from_type, "-", multi_connection_pairs$to_type)

# fix the reverse
multi_connection_pairs$connection_type[multi_connection_pairs$connection_type == "gv-lc"] <- "lc-gv"
multi_connection_pairs$connection_type[multi_connection_pairs$connection_type == "plv-lc"] <- "lc-plv"











# # try my own little grpah with only ill to see how it looks ---------------
# 
# df <- big_connection_df %>% filter(type == "gene_sharing")
# graph <- as_tbl_graph(df, directed = F)
# V(graph)$type <- contig_df$type[match(V(graph)$name, contig_df$contig_id)]
# V(graph)$color <- contig_df$color[match(V(graph)$name, contig_df$contig_id)]
# 
# 
# standart_layout <- create_layout(graph, layout = "fr")
# 
# ggraph(standart_layout) +
#   geom_edge_link(alpha = 0.5) +
#   geom_node_point(aes(colour = type))
# 
# # focus on one group
# graph <- graph %>% 
#   mutate(louvain_group = group_louvain())
# 
# target_group <- graph %>% 
#   as_tibble() %>% 
#   filter(name == "Aved_tig00303955-10-54120_vph") %>% 
#   pull(louvain_group)
# 
# 
# subgraph <- graph %>%
#   filter(louvain_group == target_group)
# 
# ggraph(subgraph, layout = "fr") +
#   geom_edge_link(alpha = 0.5) +
#   geom_node_point(aes(colour = type))
# 
# 
# # test visnetwork for nice layouts
# nodes <- data.table(id = names(V(subgraph)))
# edges <- edgelist_gene_sharing %>% filter(from %in% names(V(subgraph)), to %in% names(V(subgraph))) %>% 
#   select(from, to)
# visNetwork::visNetwork(nodes, edges)
# visNetwork::visNetwork(nodes, edges) %>% 
#   visNetwork::visIgraphLayout(layout = "layout_with_fr") %>% 
#   visNetwork::visEdges(arrows = "middle")
# 
# 
# 
# 
# 
# 
# # how many unique nodes are in all of the layers?
# for(t in unique(big_connection_df$type)){
#   print(t)
#   a <- big_connection_df %>% 
#     filter(type == t) %>% 
#     select(from, to) %>% 
#     unlist() %>% 
#     unique() %>% 
#     length()
#   print(a)
#   
# } 
# 
# 
# 
# 
# 
# 
# g_occ_ill <- graph_from_data_frame(edgelist_occurance_ill %>% filter(correlation == "positive") %>% select(from, to, spearman_ont))
# 
# V(g_occ_ill)$type <- contig_df$type[match(V(g_occ_ill)$name, contig_df$contig_id)]
# V(g_occ_ill)$color <- contig_df$color[match(V(g_occ_ill)$name, contig_df$contig_id)]
# 
# ggraph(g_occ_ill) +
#   geom_edge_link() +
#   geom_node_point()




