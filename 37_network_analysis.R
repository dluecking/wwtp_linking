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
library(googlesheets4)
library(tidygraph)
library(ggiraph)

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

# tax information
# load GV info
sheet_url <- "https://docs.google.com/spreadsheets/d/1QLNiqSt0XOS4xVPAeZAppwVjjjPIKdEE6w6f2_Qm55c/edit?gid=1228834474#gid=1228834474"
GV_info <- read_sheet(sheet_url, sheet = "Final GVs overview") %>% 
  select(shortname, personal_assessment_order)

# load tax info
lc_tax_info <- readRDS("intermediate/lc_tax/lc_tax_info_df.csv")
lc_tax_info_short <- rbindlist(lapply(list.files("intermediate/lc_tax/", pattern = ".*filtered.csv", full.names = TRUE), fread))

# add tax info to contig_df
# first all lcs:
contig_df$tax_info <- lc_tax_info_short$majority_organism[match(contig_df$contig_id, lc_tax_info_short$contig_id)]

# then mutate
contig_df <- contig_df %>% 
  mutate(
    tax_info = case_when(
      str_ends(contig_id, "vph") ~ contig_id,
      str_ends(contig_id, "plv") ~ contig_id,
      str_ends(contig_id, "lc") ~ tax_info,
      TRUE ~ "GV - to be added"
    )
  )

# then add GV
for(i in 1:nrow(contig_df)){
  if(contig_df$type[i] == "gv")
    if(length((GV_info$personal_assessment_order[GV_info == contig_df$contig_id[i]])) > 0){
      contig_df$tax_info[i] <- GV_info$personal_assessment_order[GV_info == contig_df$contig_id[i]]
    } else {
      contig_df$tax_info[i] <- "unknown"
    }
}


# load data ---------------------------------------------------------------

edgelist_crispr <- fread("intermediate/network/crispr.csv")
edgelist_integration_b <- fread("intermediate/network/integration_b.csv")
edgelist_integration_m <- fread("intermediate/network/integration_m.csv")
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing_only_interesting.csv")
edgelist_occurance_ill <- fread("intermediate/network/occurance_ill_test.csv") # ATTENTION! This one is currently set to the stricter occurance filter
edgelist_occurance_ont <- fread("intermediate/network/occurance_ont_test.csv") # ATTENTION! This one is currently set to the stricter occurance filter
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
df <- big_connection_df %>% 
  # filter(!str_detect(type, ".*occurance.*")) %>% 
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
V(g)$tax_info <- contig_df$tax_info[match(V(g)$name, contig_df$contig_id)]

E(g)$edge_type <- as.factor(E(g)$type)



g <- as_tbl_graph(g)

layout <- create_layout(g, layout = "fr")
ggiraph_plot <- ggraph(layout) +
  geom_edge_link(aes(color = type), show.legend = F) +
  geom_point_interactive(size = 3,
                         alpha = 0.8,
                         aes(x = x, 
                             y = y, 
                             color = node_type,
                             data_id = name,
                             tooltip = paste("Node ID:", name,
                                             "\nType:", node_type,
                                             "\nTaxonomy:", tax_info))) +
  scale_color_manual(values = c(
    vph = "goldenrod1",
    lc = "seagreen",
    plv = "hotpink",
    gv = "steelblue"
  )) +
  facet_edges(~edge_type) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(
    axis.text.x = element_blank(),   # Remove x-axis labels
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.ticks = element_blank(),    # Remove axis tick marks
    axis.line = element_blank()      # Remove axis lines
  )

interactive_plot <- girafe(ggobj = ggiraph_plot)
htmltools::save_html(interactive_plot, "final/subcluster_1_interactive.html")
interactive_plot

# another idea: GV-LC partners --------------------------------------------



VIRUS_HOST_df <- data.table(gv_id = as.character(),
                            host_id = as.character(),
                            gv_order = as.character(),
                            host_tax = as.character(),
                            no_of_layers = as.numeric(),
                            layers = as.character()
)



# Create VIRUS_HOST_df using data.table operations
VIRUS_HOST_df <- big_connection_df[
  # Step 1: Identify if 'from'/'to' are potential hosts or GVs
  # .SD is a special data.table symbol referring to the subset of data for a group
  , `:=`(
    is_from_host = endsWith(from, "_lc"),
    is_to_host = endsWith(to, "_lc"),
    is_from_gv = from %in% GV_info$shortname,
    is_to_gv = to %in% GV_info$shortname
  )
][
  # Step 2: Assign gv_id and host_id based on a clear GV-Host pair
  # fifelse is data.table's optimized ifelse
  , `:=`(
    gv_id = fifelse(is_from_gv & is_to_host, from,
                    fifelse(is_to_gv & is_from_host, to, NA_character_)),
    host_id = fifelse(is_from_gv & is_to_host, to,
                      fifelse(is_to_gv & is_from_host, from, NA_character_))
  )
][
  # Step 3: Filter out rows that don't represent a clear GV-Host interaction
  !is.na(gv_id) & !is.na(host_id)
][
  # Step 4: Group by gv_id and host_id and summarize types
  , .(
    layers = paste(sort(unique(type)), collapse = ", "), # Collect and sort unique types
    no_of_layers = uniqueN(type)                         # Count unique types
  ), by = .(gv_id, host_id)
][
  # Step 5: Add placeholder columns and reorder
  , `:=`(
    gv_order = NA_character_, # Will be filled later
    host_tax = NA_character_  # Will be filled later
  )
][
  # Step 6: Select and reorder final columns
  , .(gv_id, host_id, gv_order, host_tax, no_of_layers, layers)
]

VIRUS_HOST_df <- VIRUS_HOST_df %>% 
  filter(no_of_layers >= 2)

# add order information
VIRUS_HOST_df$gv_order <- GV_info$personal_assessment_order[match(VIRUS_HOST_df$gv_id, GV_info$shortname)]
VIRUS_HOST_df$host_tax <- lc_tax_info_short$majority_organism[match(VIRUS_HOST_df$host_id, lc_tax_info_short$contig_id)] 
VIRUS_HOST_df$host_tax_perc <- as.character(lc_tax_info_short$pct_reads_assigned_to_majority_taxon[match(VIRUS_HOST_df$host_id, lc_tax_info_short$contig_id)])
VIRUS_HOST_df$host_tax_perc <- str_replace(VIRUS_HOST_df$host_tax_perc, "\\.", ",")

fwrite(VIRUS_HOST_df, "tmp_V_H_pairs.csv")


# multi_connection_pairs <- big_connection_df %>%
#   filter(!str_detect(type, "negative")) %>% 
#   group_by(from, to) %>%
#   summarise(
#     unique_connection_types = n_distinct(type)
#   ) %>%
#   ungroup() %>% # Ungroup to remove grouping structure
#   filter(unique_connection_types >= 1)
# 
# 
# # now add information about type
# multi_connection_pairs$from_type <- contig_df$type[match(multi_connection_pairs$from, contig_df$contig_id)]
# multi_connection_pairs$to_type <- contig_df$type[match(multi_connection_pairs$to, contig_df$contig_id)]
# 
# 
# multi_connection_pairs$connection_type <- paste0(multi_connection_pairs$from_type, "-", multi_connection_pairs$to_type)
# 
# # fix the reverse
# multi_connection_pairs$connection_type[multi_connection_pairs$connection_type == "gv-lc"] <- "lc-gv"
# multi_connection_pairs$connection_type[multi_connection_pairs$connection_type == "plv-lc"] <- "lc-plv"
# 
# # filter for only virus-host-pairs
# virus_host_pairs <- multi_connection_pairs %>% 
#   filter(connection_type == "lc-gv")
# 
# 
# GVs_with_host <- unique(c(virus_host_pairs$from, virus_host_pairs$to))
# GVs_with_host <- GVs_with_host[!grepl("lc$", GVs_with_host)]




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




