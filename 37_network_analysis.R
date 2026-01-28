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
edgelist_occurance_ill <- fread("intermediate/network/occurance_ill.csv") 
edgelist_occurance_ont <- fread("intermediate/network/occurance_ont.csv")
edgelist_non_crispr <- fread("intermediate/network/non_CRISPR.csv")


# load crispr cas data ----------------------------------------------------

crispr_df <- fread("intermediate/CRISPR/cassette/HMM2019_cassettes.csv") %>% 
  filter(bitscore >= 100)
crispr_df$contig_id <- str_remove(crispr_df$V1, "\\_\\d+\\_ID.*$")


# explore -----------------------------------------------------------------
# create binary yes/no big list with connection type for each edge
big_connection_df <- rbind(
  edgelist_crispr %>% select(from, to, value = "crispr") %>% 
    mutate(type = "crispr"),
  edgelist_gene_sharing %>% select(from, to, value = "gene_sharing") %>% 
    mutate(type = "gene_sharing"),
  edgelist_integration_b %>% select(from, to, value = "integration_b") %>% 
    mutate(type = "integration_boundary"),
  edgelist_integration_m %>% select(from, to, value = "integration_m") %>% 
    mutate(type = "integration_middle"),
  edgelist_non_crispr %>% select(from, to, value = "non_CRISPR") %>% 
    mutate(type = "non_crispr"),
  edgelist_occurance_ill %>% filter(correlation == "positive") %>% select(from, to, value = "absolut_spearman") %>% 
    mutate(type = "occurance_ill_positive"),
  edgelist_occurance_ont %>% filter(correlation == "positive") %>% select(from, to, value = "absolut_spearman") %>% 
    mutate(type = "occurance_ont_positive"),
  edgelist_occurance_ill %>% filter(correlation == "negative") %>% select(from, to, value = "absolut_spearman") %>% 
    mutate(type = "occurance_ill_negative"),
  edgelist_occurance_ont %>% filter(correlation == "negative") %>% select(from, to, value = "absolut_spearman") %>% 
    mutate(type = "occurance_ont_negative")
)

rm(edgelist_crispr, edgelist_gene_sharing, edgelist_integration_b, edgelist_integration_m, edgelist_non_crispr, edgelist_occurance_ill, edgelist_occurance_ont)



# visualize a subcluster surrounding a specific node ----------------------

CONTIG_OF_INTEREST <- "AalE_tig00021708-10-192480_vph" # thats the good one
# CONTIG_OF_INTEREST <- "Bjer_2_3"
SAVE_PLOT <- FALSE

plvs <- str_remove(list.files("intermediate/contigs/plv"), "\\.fna")
vphs <- str_remove(list.files("intermediate/contigs/vph"), "\\.fna")
gvs <- GV_info$shortname

list_of_sequences_to_print <- vphs

# pre calc for each one the same:
# I had this code chunk, I cant remember why
# df <- big_connection_df %>% 
#   filter(!str_detect(type, ".*occurance.*")) %>% 
#   filter(type != "integration_middle")

graph <- as_tbl_graph(big_connection_df, directed = F)
graph <- graph %>% 
  mutate(louvain_group = group_louvain())

for(contig in list_of_sequences_to_print){
  CONTIG_OF_INTEREST <- contig
  
  # get the louvain group of the current contig
  target_group <- graph %>% 
    as_tibble() %>% 
    filter(name == CONTIG_OF_INTEREST) %>%
    pull(louvain_group)
  
  # if this is not in any subcluster, we can skip it
  if(identical(target_group, integer(0))){
    next
  }
  
  # filter the graph to contain only this louvain group
  subgraph <- graph %>%
    filter(louvain_group == target_group)
  
  # construct df that contains info on all contigs in this subgraph
  small_node_df <- data.table(id = names(V(subgraph)),
                              contig_type = "")
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
  small_edge_df <- big_connection_df %>%
    filter(from %in% small_node_df$id & to %in% small_node_df$id)
  
  g <- graph_from_data_frame(small_edge_df, directed = F, small_node_df)
  V(g)$node_type <- contig_df$type[match(V(g)$name, contig_df$contig_id)]
  V(g)$tax_info <- contig_df$tax_info[match(V(g)$name, contig_df$contig_id)]
  
  E(g)$edge_type <- as.factor(E(g)$type)
  
  
  g <- as_tbl_graph(g)
  
  # CHANGE HERE!
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
  interactive_plot
  
  
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

g <- graph_from_data_frame(big_connection_df)

g <- graph_from_data_frame(d = big_connection_df, 
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

# 2. Plot with ggplot2
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




# does contig length correlate with K (number of connections?) ------------

deg <- degree(g, mode = "all")




# centrality --------------------------------------------------------------

# we need to do this for each layer seperately and then one combined:
for(layer in c(unique(big_connection_df$type), "all")){
  print(layer)
  
  if(layer == "all"){
    g <- graph_from_data_frame(big_connection_df)
  }else{
    g <- graph_from_data_frame(big_connection_df %>% filter(type == layer))
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

# THIS PART IS EXPLORATORY ##########
SUBSAMPLE_SIZE <- 0.1
edgelist_gene_sharing <- fread("intermediate/network/gene_sharing.csv")
edgelist_gene_sharing <- edgelist_gene_sharing %>% 
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

big_connection_df <- edgelist_gene_sharing %>% select(from, to) %>% 
  mutate(type = "gene_sharing") %>% 
  sample_n(size = nrow(edgelist_gene_sharing) * SUBSAMPLE_SIZE)



g <- graph_from_data_frame(big_connection_df %>% filter(type == "gene_sharing"))

deg <- degree(g, mode = "all")
btw <- betweenness(g, directed = FALSE, normalized = TRUE)

deg_df <- enframe(deg, name = "contig_id", value = "degree")
btw_df <- enframe(btw, name = "contig_id", value = "betweeness")

network_df <- left_join(deg_df, btw_df)
network_df$contig_type <- contig_df$type[match(network_df$contig_id, contig_df$contig_id)]

# add the information: is the lc connected to a gv or not?
# for(i in 1:nrow(network_df)){
#   if(i %% 100 == 0){
#     print(i)
#   }
#   current_lc <- network_df$contig_id[i]
#   
#   # this only applies to LCs, skip if non-lc
#   if(!str_ends(current_lc, "\\_lc$")){
#     next
#   }
#   
#   # else check if we are connected to a GV though gene sharing
#   tmp_df <- edgelist_gene_sharing %>% 
#     filter(from == current_lc | to == current_lc)
#   
#   if(any(tmp_df$from_type == "gv" | tmp_df$to_type == "gv")){
#     network_df$contig_type[i] <- "lc_gv_connected"
#   }
#   
# }

# First find all LC contigs that are connected to a GV
lc_gv_connected <- edgelist_gene_sharing %>%
  filter(from_type == "gv" | to_type == "gv") %>%   # only GV edges
  mutate(lc = case_when(
    str_ends(from, "_lc") ~ from,
    str_ends(to, "_lc")   ~ to,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(lc)) %>%
  distinct(lc) %>%
  pull(lc)

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
    map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15,) +
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
    map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15,) +
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


# gene sharing centrality -------------------------------------------------
# the idea is: lcs which are connected to GVs are on average connected more to other entities, than LCs that are not connected to GVs
# this goes to SUPP, but the idea is built in the figure above
# 
# edgelist_gene_sharing <- fread("intermediate/network/gene_sharing_only_interesting.csv")
# 
# all_strings <- c(edgelist_gene_sharing$from, edgelist_gene_sharing$to)
# all_lcs <- unique(all_strings[str_ends(all_strings, "lc")])
# 
# lc_df <- data.table(contig = all_lcs,
#                     connections = 0,
#                     is_GV_connected = FALSE)
# 
# for(i in 1:nrow(lc_df)){
#   # which lc are we looking at?
#   current_LC <- lc_df$contig[i]
#   
#   # filter to only retain tmp with this lc
#   tmp_df <- edgelist_gene_sharing %>% 
#     filter(from == current_LC | to == current_LC)
#   
#   # is there a gv connected to the lc?
#   if(any(tmp_df$from_type == "gv" | tmp_df$to_type == "gv")){
#     lc_df$is_GV_connected[i] <- TRUE
#   }
#   
#   # how many unique connections do we count?
#   tmp_df <- tmp_df %>%
#     rowwise() %>%
#     mutate(connection = paste(sort(c(from, to)), collapse = "-"))
#   
#   lc_df$connections[i]  <- n_distinct(tmp_df$connection)
# }
# 
# p <- ggplot(lc_df, aes(x = is_GV_connected, y = connections, fill = is_GV_connected)) +
#   geom_signif(
#     comparisons = list(c("FALSE", "TRUE")),
#     map_signif_level = TRUE, textsize = 3, step_increase = 0.1, margin_top = 0.15) +
#   geom_boxplot() +
#   theme_cowplot() +
#   theme(legend.position = "None") +
#   ylim(c(0,NA)) +
#   scale_y_log10() +
#   scale_x_discrete(labels = c("not-connected", "GV-connected")) +
#   scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "steelblue")) +
#   ggtitle(label = "Gene Sharing Comparison",
#           subtitle = "Number of connections of gv-connected vs non-connected LCs") +
#   ylab("# of gene sharing connections") +
#   xlab("")
# 
# ggsave(plot = p, file = "final/centrality_plots/gv_vs_non-gv-connected_LCs.svg", height = 6, width = 6)
# ggsave(plot = p, file = "final/centrality_plots/gv_vs_non-gv-connected_LCs.pdf", height = 6, width = 6)
# ggsave(plot = p, file = "final/centrality_plots/gv_vs_non-gv-connected_LCs.png", height = 6, width = 6)




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



# PLVs and VPHs -----------------------------------------------------------

plv_df <- big_connection_df %>% 
  filter(str_detect(from, "vph") | str_detect(to, "vph"))



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




