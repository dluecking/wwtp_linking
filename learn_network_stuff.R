# Author: dlu @ veelab
# Version: 2025-06-17

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyr)
library(patchwork)
library(cowplot)
library(igraph)
library(igraphdata)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load data ---------------------------------------------------------------

edgelist_crispr <- fread("intermediate/network/crispr.csv")
edgelist_integration_b <- fread("intermediate/network/integration_b.csv")
edgelist_integration_m <- fread("intermediate/network/integration_m.csv")
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing_only_interesting.csv")

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



# exploring ---------------------------------------------------------------
# create network plot from the integration_b layer 
g_int_b <- graph_from_data_frame(edgelist_integration_b, directed = F)
# this assigns a type to each node in the network
V(g_int_b)$type <- contig_df$type[match(V(g_int_b)$name, contig_df$contig_id)]

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

# add to graph
V(g_int_b)$color <- contig_df$color[match(V(g_int_b)$name, contig_df$contig_id)]

g_int_b
plot(g_int_b,edge.arrow.size=1, vertex.size=7, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.5, vertex.label.dist=1, edge.curved=0.1,layout=layout_with_fr)

# CENTRALITY (we test this for gene sharing)
g_gene_sharing <- graph_from_data_frame(edgelist_gene_sharing, directed = T)


centr_degree(g_gene_sharing, mode = "all")

# this assigns a type to each node in the network
V(g_gene_sharing)$type <- contig_df$type[match(V(g_gene_sharing)$name, contig_df$contig_id)]
# add to graph
V(g_gene_sharing)$color <- contig_df$color[match(V(g_gene_sharing)$name, contig_df$contig_id)]

plot(g_gene_sharing,edge.arrow.size=.01, vertex.size=3, 
     vertex.frame.color="black", vertex.label = NA, 
     vertex.label.cex=.5, vertex.label.dist=1, edge.curved=0.1,layout=layout_with_fr)


E(g_gene_sharing)

# NOTES
cluster_walktrap(g_gene_sharing)
cluster_louvain(g_gene_sharing)




