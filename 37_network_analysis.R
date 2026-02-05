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
library(patchwork)
library(tibble)
library(ggsignif)
library(doParallel)
library(foreach)

# setup -------------------------------------------------------------------

# cores
cores_to_use <- 14
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

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
  select(shortname, personal_assessment_order, public_ID) %>% 
  drop_na(shortname)

# load other info
sheet_url <- "https://docs.google.com/spreadsheets/d/1CnqcfhOfS0rVBU6mvyCKrm50BmxJjVvif47f9vkQZgU/edit?gid=836958150#gid=836958150"
vph_info <- read_sheet(sheet_url, sheet = "Table S3") %>% 
  select(contig_ID, public_ID)
plv_info <- read_sheet(sheet_url, sheet = "Table S4") %>% 
  select(contig_ID, public_ID)
vph_plv_combined_info <- rbind(vph_info, plv_info)

# load tax info
lc_tax_info <- readRDS("intermediate/lc_tax/lc_tax_info_df.csv")
lc_tax_info_short <- rbindlist(lapply(list.files("intermediate/lc_tax/", pattern = ".*filtered.csv", full.names = TRUE), fread))

# add tax info to contig_df
# first all lcs (which get a (76%) tag, which indicates the percentage of genes mapped to this tax)
contig_df$tax_info <- paste0(lc_tax_info_short$majority_organism[match(contig_df$contig_id, lc_tax_info_short$contig_id)],
                             " (", round(lc_tax_info_short$pct_reads_assigned_to_majority_taxon[match(contig_df$contig_id, lc_tax_info_short$contig_id)], 2), " %)")

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
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing.csv")
edgelist_occurance_ill <- fread("intermediate/network/occurance_ill_3.csv") 
edgelist_occurance_ont <- fread("intermediate/network/occurance_ont_3.csv")
edgelist_non_crispr <- fread("intermediate/network/non_CRISPR.csv")

# crispr data
crispr_df <- fread("intermediate/CRISPR/cassette/HMM2019_cassettes.csv") %>% 
  filter(bitscore >= 100)
crispr_df$contig_id <- str_remove(crispr_df$V1, "\\_\\d+\\_ID.*$")


# quick explore to get to overlap of ill and ont --------------------------
# add ill edge id and clean df
edgelist_occurance_ill <- edgelist_occurance_ill %>%
  mutate(
    node1 = pmin(from, to),
    node2 = pmax(from, to),
    edge_id = paste0(node1, "--", node2)
  ) %>%
  select(edge_id, node1, node2, 
         spearman_ill = spearman_ont, # this is a mistake in another script
         correlation_ill = correlation) %>%
  distinct(edge_id, .keep_all = TRUE)

# same for ont
edgelist_occurance_ont <- edgelist_occurance_ont %>%
  mutate(
    node1 = pmin(from, to),
    node2 = pmax(from, to),
    edge_id = paste0(node1, "--", node2)
  ) %>%
  select(edge_id, node1, node2, 
         spearman_ont = spearman_ont, 
         correlation_ont = correlation) %>%
  distinct(edge_id, .keep_all = TRUE)

# inner join to get only shared edges
occurance_edges <- edgelist_occurance_ill %>%
  inner_join(edgelist_occurance_ont, by = c("edge_id", "node1", "node2")) %>%
  mutate(
    agree = (correlation_ill == correlation_ont)
  ) %>%
  select(edge_id, node1, node2, 
         spearman_ill, correlation_ill, 
         spearman_ont, correlation_ont, 
         agree)

# remove where we are in disagreement (3 edges out of 1.5 million)
occurance_edges <- occurance_edges %>% filter(agree == "TRUE")
# create this weird paspte, actually might not be necessary
occurance_edges$value <- paste0(round(occurance_edges$spearman_ill, 2), "_ill--", round(occurance_edges$spearman_ont, 2), "_ont")
# rename
names(occurance_edges) <- c("edge_id", "from", "to", "spearman_ill", "correlation_ill", "spearman_ont", "correlation_ont", "agree", "value") 

a <- occurance_edges %>%
  filter(!grepl("_lc$", from) | !grepl("_lc$", to))

# prepare single df -------------------------------------------------------
# create large combined df, edgelist with weighted values
big_connection_df <- rbind(
  edgelist_crispr %>% select(from, to, value = "crispr") %>% 
    mutate(type = "crispr"),
  edgelist_gene_sharing %>% 
    # filter(gene_sharing >= 0.05) %>% # We do this in the other script that prepares gene sharing
    select(from, to, value = "gene_sharing") %>% 
    mutate(type = "gene_sharing"),
  edgelist_integration_b %>% select(from, to, value = "integration_b") %>% 
    mutate(type = "integration_boundary"),
  edgelist_integration_m %>% select(from, to, value = "integration_m") %>% 
    mutate(type = "integration_middle"),
  edgelist_non_crispr %>% select(from, to, value = "non_CRISPR") %>% 
    mutate(type = "non_crispr"),
  occurance_edges %>%
    filter(correlation_ont == "positive") %>% 
    select(from, to, value = value) %>% 
    mutate(type = "occurance_positive"),
  occurance_edges %>%
    filter(correlation_ont == "negative") %>% 
    select(from, to, value = value) %>% 
    mutate(type = "occurance_negative")
)

rm(edgelist_crispr, edgelist_gene_sharing, edgelist_integration_b, edgelist_integration_m, edgelist_non_crispr, edgelist_occurance_ill, edgelist_occurance_ont)
rm(occurance_edges)

# create unique edge id
big_connection_df_filtered <- big_connection_df %>%
  mutate(edge_id = paste0(pmin(from, to), "--", pmax(from, to)))

# keep only one row per edge_id + type (removes duplicate connections but we lose info on weightedness)
big_connection_df_filtered <- big_connection_df_filtered %>% 
  dplyr::distinct(edge_id, type, .keep_all = TRUE)

rm(big_connection_df)


# test weird community structure ------------------------------------------

res <- list()

for(t in c(unique(big_connection_df_filtered$type), "none")){
  if(t != "none"){
    graph <- as_tbl_graph(big_connection_df_filtered %>% filter(type != t), directed = FALSE)
  } else {
    graph <- as_tbl_graph(big_connection_df_filtered, directed = FALSE)
  }
  graph <- graph %>%
    mutate(louvain_group = group_louvain())
  a <- as.data.table(table(V(graph)$louvain_group))
  setnames(a, c("V1", "N"), c("group", "count"))
  a[, group := as.integer(group)]
  a[, removed_type := t]
  res[[t]] <- a
}
dt <- rbindlist(res)

ecdf_plot <- ggplot(dt, aes(x = count, colour = removed_type)) +
  stat_ecdf() +
  scale_x_log10() +
  theme_minimal() +
  labs(
    x = "Community size (log scale)",
    y = "ECDF",
    title = "Distribution of Louvain community sizes"
  )
ggsave(plot = ecdf_plot, filename = "testing/louvain_plots/ecdf.png")

multi_plot <- ggplot(dt, aes(x = group, y = count)) +
  geom_col() +
  facet_wrap(~ removed_type, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Louvain group",
    y = "Number of contigs",
    title = "Louvain community size distributions",
    subtitle = "Panel = this has been removed"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(plot = multi_plot, filename = "testing/louvain_plots/multipanel.png")


# test resolutions of network:
for(res in c(0.5, 1.0, 1.5, 2.0, 3.0)) {
  communities <- cluster_louvain(graph, resolution = res)
  sizes <- table(membership(communities))
  
  cat("\nResolution:", res, "\n")
  cat("  N communities:", length(sizes), "\n")
  cat("  Largest:", max(sizes), "nodes\n")
  cat("  Communities 10-500 nodes:", sum(sizes >= 10 & sizes <= 500), "\n")
  cat("  Singletons:", sum(sizes == 1), "\n")
}

# visualize a subcluster surrounding a specific node ----------------------

# CONTIG_OF_INTEREST <- "AalE_tig00021708-10-192480_vph" # thats the good one
CONTIG_OF_INTEREST <- "Bjer_2_3"
SAVE_PLOT <- TRUE

# this is only for batch creation of plots:
plvs <- str_remove(list.files("intermediate/contigs/plv"), "\\.fna")
vphs <- str_remove(list.files("intermediate/contigs/vph"), "\\.fna")
gvs <- GV_info$shortname

list_of_sequences_to_print <- plvs


# this is creating the graph and groups
graph <- as_tbl_graph(big_connection_df_filtered, directed = F)
graph <- graph %>% 
  mutate(louvain_group = group_louvain(resolution = 3))


for(contig in list_of_sequences_to_print){
  CONTIG_OF_INTEREST <- contig
  
  # get the louvain group of the current contig
  target_group <- graph %>% 
    as_tibble() %>% 
    filter(name == CONTIG_OF_INTEREST) %>%
    pull(louvain_group)
  
  # if this is not in any subcluster, we can skip it
  if(identical(target_group, integer(0))){
    cat("This one is NOT in a group:", CONTIG_OF_INTEREST, "\n")
    next
  }
  
  # filter the graph to contain only this louvain group
  subgraph <- graph %>%
    filter(louvain_group == target_group)
  
  # construct df that contains info on all contigs in this subgraph
  small_node_df <- data.table(id = names(V(subgraph)),
                              contig_type = "")
  
  
  # if this subgraph is too large (1k for now) then skip
  if(nrow(small_node_df) >= 500){
    cat("This one was too big:", CONTIG_OF_INTEREST, "\n")
    cat("With this many nodes:", nrow(small_node_df), "\n")
    next
  }
    
    
    
  small_node_df$contig_type <- contig_df$type[match(small_node_df$id, contig_df$contig_id)]
  small_node_df$public_ID <- small_node_df$id
  small_node_df$public_ID[small_node_df$contig_type == "gv"] <- 
    GV_info$public_ID[match(small_node_df$id[small_node_df$contig_type == "gv"], GV_info$shortname)]
  
  # CRISP info needs to be added
  small_node_df$cas_genes <- ""
  for(i in 1:nrow(small_node_df)){
    # virophages get an NA
    if(small_node_df$contig_type[i] == "vph"){
      small_node_df$cas_genes[i] <- "Not applicable"
      next
    }
    
    # GVs and LCs get a concatenated unique list of cas genes found
    tmp_df<- crispr_df %>% 
      filter(contig_id == small_node_df$id[i])
    
    if(nrow(tmp_df) >= 1){
      cas_genes <- tmp_df %>% 
        select(annotation) %>% 
        unlist() %>% 
        unique() %>% 
        paste(collapse = ", ")
    }else{
      cas_genes <- "None detected."
    }
    small_node_df$cas_genes[i] <- cas_genes
  }
  
  # then subset the big df, so we retain only edges of nodes that are part of the 
  # group
  small_edge_df <- big_connection_df_filtered %>%
    filter(from %in% small_node_df$id & to %in% small_node_df$id)
  
  g <- graph_from_data_frame(small_edge_df, directed = F, small_node_df)
  V(g)$node_type <- contig_df$type[match(V(g)$name, contig_df$contig_id)]
  V(g)$tax_info <- contig_df$tax_info[match(V(g)$name, contig_df$contig_id)]
  
  E(g)$edge_type <- as.factor(E(g)$type)
  
  
  g <- as_tbl_graph(g)
  
  # highlight the current contig int the plot
  g <- g %>%
    activate(nodes) %>%
    mutate(is_focal = (name == CONTIG_OF_INTEREST))
  
  # create layout
  layout <- create_layout(g, layout = "fr")
  # plot
  ggiraph_plot <- ggraph(layout) +
    geom_edge_link(aes(color = type), show.legend = F) +
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
                                               "\nCas genes: ", cas_genes),
                               color = if_else(is_focal, "black", "white"),   # black border for focal node, none for others
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
    guides(fill = guide_legend(
      override.aes = list(size = 4)  # increase legend point size
    ))+
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
        "occurance_ill_positive" = "Ill. positive", 
        "occurance_ont_positive" = "ONT positive", 
        "occurance_ill_negative" = "Ill. negative", 
        "occurance_ont_negative" = "ONT negative"
      )
    )
    ) +
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
  # interactive_plot
  
  
  if(SAVE_PLOT){
    type_of_contig <- str_remove(CONTIG_OF_INTEREST, "^.*\\_")
    
    if(type_of_contig == "plv" || type_of_contig == "vph"){
      public_ID <- vph_plv_combined_info$public_ID[vph_plv_combined_info$contig_ID == CONTIG_OF_INTEREST]
      FILE_BASE <- paste0("final/plv_vph_subclusters/", public_ID) 
    }else{
      public_ID <- GV_info$public_ID[GV_info$shortname == CONTIG_OF_INTEREST]
      FILE_BASE <- paste0("final/gv_subclusters/", public_ID)
    }
    
    # in case you want to save
    # for ggsave this is a bit more complicated, since we need to adjust the width and height:
    FACETS <- length(unique(small_edge_df$type))
    HEIGHT_FACTOR <- ceiling(FACETS / 3)
    
    ggsave(ggiraph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".svg"), 
           width = 4.5, height = 1.8*HEIGHT_FACTOR)
    ggsave(ggiraph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".pdf"), 
           width = 4.5, height = 1.8*HEIGHT_FACTOR)
    ggsave(ggiraph_plot + theme(legend.position = "none"), 
           file = paste0(FILE_BASE, ".png"), 
           width = 4.5, height = 1.8*HEIGHT_FACTOR)
    htmltools::save_html(interactive_plot, file = paste0(FILE_BASE, "_interactive.html"))
  }
}



# check scale free and other stats ----------------------------------------

g <- graph_from_data_frame(d = big_connection_df_filtered, 
                           directed = FALSE, 
                           vertices = contig_df)

# 1. Generate the data
d <- degree(g)

# Create a frequency table and calculate the cumulative probability
degree_counts <- as.data.frame(table(d))
names(degree_counts) <- c("k", "Freq")
degree_counts$k <- as.numeric(as.character(degree_counts$k))

# Sort and calculate cumulative distribution P(X >= k)
degree_counts <- degree_counts[order(degree_counts$k, decreasing = TRUE), ]
degree_counts$cumulative_prob <- cumsum(degree_counts$Freq) / sum(degree_counts$Freq)

# plot
ggplot(degree_counts, aes(x = k, y = cumulative_prob)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  scale_x_log10() + 
  scale_y_log10() +
  annotation_logticks() + # Adds the visual "rug" lines for log scales
  labs(
    title = "Log-Log Degree Distribution",
    subtitle = "Check for linearity to identify scale-free properties",
    x = "Degree (k)",
    y = "Cumulative Probability P(k)"
  ) +
  theme_minimal()



# check with poweRlaw -----------------------------------------------------

library(poweRlaw)

# Create a distribution object
d_positive <- d[d > 0]
m_pl <- displ$new(d_positive)

# Estimate the minimum degree (xmin) where the power law starts
est <- estimate_xmin(m_pl)
m_pl$setXmin(est)

# Compare with a log-normal distribution
m_ln <- dislnorm$new(d_positive)
m_ln$setXmin(m_pl$getXmin()) # Compare on the same tail
m_ln$setPars(estimate_pars(m_ln))

comp <- compare_distributions(m_pl, m_ln)
comp$test_statistic # Positive = Power law is better; Negative = Log-normal is better
comp$p_two_sided

rm(d, d_positive, g)


# centrality --------------------------------------------------------------

g <- graph_from_data_frame(d = big_connection_df_filtered, 
                           directed = FALSE, 
                           vertices = contig_df)

# if you want a single layer
layer <- "all"

# we need to do this for each layer seperately and then one combined:
for(layer in c(unique(big_connection_df_filtered$type), "all")){
  print(layer)
  
  if(layer == "all"){
    g <- graph_from_data_frame(big_connection_df_filtered)
  }else{
    g <- graph_from_data_frame(big_connection_df_filtered %>% filter(type == layer))
  }
  
  deg <- degree(g, mode = "all")
  btw <- betweenness(g, directed = FALSE, normalized = TRUE)
  
  deg_df <- enframe(deg, name = "contig_id", value = "degree")
  btw_df <- enframe(btw, name = "contig_id", value = "betweeness")
  
  network_df <- left_join(deg_df, btw_df)
  network_df$contig_type <- contig_df$type[match(network_df$contig_id, contig_df$contig_id)]
  
  # set order of entities
  network_df$contig_type <- factor(network_df$contig_type,
                                   levels = c("lc", "gv", "plv", "vph"))
  
  deg_plot <- ggplot(network_df, aes(x = contig_type, y = degree, fill = contig_type)) +
    geom_boxplot() +
    geom_signif(
      comparisons = list(c("lc", "vph"), c("lc", "plv"), c("lc", "gv")),
      map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15,) +
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
  deg_plot  
  
  
  
  bet_plot <- ggplot(network_df, aes(x = contig_type, y = betweeness, fill = contig_type)) +
    geom_boxplot() +
    geom_signif(
      comparisons = list(c("lc", "vph"), c("lc", "plv"), c("lc", "gv")),
      map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15,) +
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
  
  p <- deg_plot + bet_plot
  ggsave(plot = p, file = paste0("final/centrality_plots/centrality_", layer, "_plot.png"), width = 7, height = 4)
  ggsave(plot = p, file = paste0("final/centrality_plots/centrality_", layer, "_plot.pdf"), width = 7, height = 4)
  ggsave(plot = p, file = paste0("final/centrality_plots/centrality_", layer, "_plot.svg"), width = 7, height = 4)
}


# centrality Figure for maintext ------------------------------------------

# SUBSAMPLE_SIZE <- 0.1
# edgelist_gene_sharing <- fread("intermediate/network/gene_sharing.csv")
# edgelist_gene_sharing <- edgelist_gene_sharing %>% 
#   mutate(from_type = case_when(
#     str_ends(from, "lc")  ~ "lc",
#     str_ends(from, "vph") ~ "vph",
#     str_ends(from, "plv") ~ "plv",
#     TRUE ~ "gv"
#   )) %>% 
#   mutate(to_type = case_when(
#     str_ends(to, "lc")  ~ "lc",
#     str_ends(to, "vph") ~ "vph",
#     str_ends(to, "plv") ~ "plv",
#     TRUE ~ "gv"
#   ))
# 
# big_connection_df <- edgelist_gene_sharing %>% select(from, to) %>% 
#   mutate(type = "gene_sharing") %>% 
#   sample_n(size = nrow(edgelist_gene_sharing) * SUBSAMPLE_SIZE)


# add contig type info
big_connection_df_filtered <- big_connection_df_filtered %>% 
  mutate(from_type = case_when(
    str_ends(from, "lc")  ~ "lc",
    str_ends(from, "vph") ~ "vph",
    str_ends(from, "plv") ~ "plv",
    TRUE ~ "gv"
  )) %>% 
  mutate(to_type = case_when(
    str_ends(to, "lc")  ~ "lc",
    str_ends(to, "vph") ~ "vph",
    str_ends(to, "plv") ~ "plv",
    TRUE ~ "gv"
  ))

# find all LC contigs that are connected to a GV
lc_gv_connected <- big_connection_df_filtered %>%
  filter(from_type == "gv" | to_type == "gv") %>%   # only GV edges
  mutate(lc = case_when(
    str_ends(from, "_lc") ~ from,
    str_ends(to, "_lc")   ~ to,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(lc)) %>%
  distinct(lc) %>%
  pull(lc)

g <- graph_from_data_frame(big_connection_df_filtered)

deg <- degree(g, mode = "all")
btw <- betweenness(g, directed = FALSE, normalized = TRUE, cutoff = 5)


# new estimate of betweenness using multiple cores and estimate only
# Split vertices into chunks for parallel processing
n_vertices <- vcount(g)
vertex_chunks <- split(1:n_vertices, cut(1:n_vertices, breaks = cores_to_use, labels = FALSE))

# Parallel computation
btw_list <- foreach(chunk = vertex_chunks, 
                    .packages = "igraph",
                    .combine = c) %dopar% {
                      estimate_betweenness(g, v = chunk, directed = FALSE, cutoff = 5)
                    }

# btw_list now contains betweenness values for all vertices
btw <- btw_list


deg_df <- enframe(deg, name = "contig_id", value = "degree")
btw_df <- enframe(btw, name = "contig_id", value = "betweeness")

network_df <- left_join(deg_df, btw_df)
network_df$contig_type <- contig_df$type[match(network_df$contig_id, contig_df$contig_id)]

# Update network_df$contig_type if contig_id is in that set
network_df <- network_df %>%
  mutate(contig_type = case_when(
    str_ends(contig_id, "_lc") & contig_id %in% lc_gv_connected ~ "lc_gv_connected",
    TRUE ~ contig_type
  ))

# set order of entities
network_df$contig_type <- factor(network_df$contig_type,
                                 levels = c("lc", "lc_gv_connected", "gv", "plv", "vph"))

deg_plot <- ggplot(network_df, aes(x = contig_type, y = degree, fill = contig_type)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("lc", "lc_gv_connected"), c("lc", "vph"), c("lc", "plv"), c("lc", "gv")),
    map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
  scale_y_log10() +
  labs(
    title = "Degree Distribution",
    subtitle = paste0("Layer: Gene Sharing"), 
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
deg_plot  


bet_plot <- ggplot(network_df, aes(x = contig_type, y = betweeness, fill = contig_type)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("lc", "lc_gv_connected"), c("lc", "gv")),
    map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
  scale_y_log10() +
  labs(
    title = "Betweenness Distribution",
    subtitle = paste0("Layer: Gene Sharing"),
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

p <- deg_plot + bet_plot
ggsave(plot = p, file = "final/centrality_plots/gene_sharing_lc_vs_lc-gv-connected_vs_others.png", width = 7, height = 4)
ggsave(plot = p, file = "final/centrality_plots/gene_sharing_lc_vs_lc-gv-connected_vs_others.pdf", width = 7, height = 4)
ggsave(plot = p, file = "final/centrality_plots/gene_sharing_lc_vs_lc-gv-connected_vs_others.svg", width = 7, height = 4)





# OLD NOT IMPORTANT
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









