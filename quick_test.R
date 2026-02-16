# After running 37b
library(dplyr)
library(data.table)
library(stringr)
library(igraph)

# setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")
setwd("/run/user/1000/gvfs/sftp:host=login01.lisc.univie.ac.at,user=luecking/lisc/home/user/luecking/luecking_scratch/projects/wwtp_linking")

big_connection_df_filtered <- fread("intermediate/network/network_analysis/big_connection_df_filtered.csv")
contig_df <- fread("intermediate/network/network_analysis/contig_df.csv")

graph <- graph_from_data_frame(
  d = big_connection_df_filtered, 
  directed = FALSE, 
  vertices = contig_df
)

communities <- cluster_louvain(graph, weights = E(graph)$weight)

membership_df <- data.frame(
  contig_id = V(graph)$name,
  louvain_group = membership(communities),
  stringsAsFactors = FALSE
)

clusters_with_mobile <- membership_df %>%
  left_join(contig_df, by = "contig_id") %>%
  group_by(louvain_group) %>%
  summarize(
    total = n(),
    n_gv = sum(type == "gv"),
    n_vph = sum(type == "vph"),
    n_plv = sum(type == "plv"),
    n_lc = sum(type == "lc")
  ) %>%
  filter(n_gv > 0 | n_vph > 0 | n_plv > 0)

nice_clusters <- clusters_with_mobile %>%
  filter(total >= 5 & total <= 150)

cat("Nice clusters with mobile elements:", nrow(nice_clusters), "\n")
print(nice_clusters)



# After loading big_connection_df_filtered
mobile_edge_counts <- big_connection_df_filtered %>%
  mutate(
    has_mobile = str_ends(from, "gv|vph|plv") | str_ends(to, "gv|vph|plv")
  ) %>%
  group_by(type, has_mobile) %>%
  summarize(n = n())

print(mobile_edge_counts)

# How many edges involve each mobile element?
mobile_nodes <- contig_df %>% filter(type %in% c("gv", "vph", "plv"))

for(mobile_id in mobile_nodes$contig_id[1:20]) {  # Check first 20
  n_edges <- sum(big_connection_df_filtered$from == mobile_id | 
                   big_connection_df_filtered$to == mobile_id)
  cat(mobile_id, ":", n_edges, "edges\n")
}


# After running 37b
library(igraph)
library(data.table)

g <- graph_from_data_frame(big_connection_df_filtered, directed = FALSE, vertices = contig_df)

# Component analysis
components <- components(g)

comp_sizes <- data.frame(
  comp_id = 1:components$no,
  size = as.numeric(table(components$membership))
) %>%
  arrange(desc(size))

cat("=== COMPONENT ANALYSIS ===\n")
cat("Total components:", components$no, "\n")
cat("Largest component:", comp_sizes$size[1], "nodes\n")
cat("Second largest:", comp_sizes$size[2], "nodes\n")
cat("Singletons:", sum(comp_sizes$size == 1), "\n\n")

cat("Percentage in largest component:", 
    round(comp_sizes$size[1] / vcount(g) * 100, 1), "%\n")
