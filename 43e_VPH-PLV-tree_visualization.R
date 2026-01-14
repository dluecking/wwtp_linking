# Author: dlu @ veelab
# Version: 2025-07-17

# Packages
library(dplyr)
library(seqinr)
library(data.table)
library(stringr)
library(ggtree)
library(ggtreeExtra)
library(googlesheets4)
library(ggnewscale)
library(tidyr)
library(ggplot2)
library(ape)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load auxiliary data -----------------------------------------------------

reference_df <- fread("helperfiles/Virophage_Table_S1_Roux2021.csv") %>% 
  filter(`Predicted completeness` == "Complete") %>% 
  filter(str_detect(`Source dataset`, pattern = "Hackl_2021|IMG_VR|NCBI Nucleotide database"))

reference_df <- reference_df %>% 
  select(`Genome Abbreviated ID`, `Genome full ID`, `Blast-based affiliation`, `Genome length`) %>% 
  mutate(
    short_id = `Genome Abbreviated ID`,
    long_id = `Genome full ID`,
    split = str_split(`Blast-based affiliation`, ";"),
    class = "Virophaviricetes",
    order = sapply(split, function(x) if (length(x) >= 2) x[2] else NA),
    family  = sapply(split, function(x) if (length(x) >= 3) x[3] else NA),
    length = `Genome length`
  ) %>% 
  select(short_id, long_id, class, order, family,length)


vph_df <- read_sheet("https://docs.google.com/spreadsheets/d/113hSsqFV73bfdHTs5WoFFQU1DvnAcV5kPYmtanhROsY/edit?gid=0#gid=0", sheet = "virophages")
vph_df <- vph_df %>% 
  mutate(
    split = str_split(taxonomy, ";"),
    class = "Virophaviricetes",
    order = sapply(split, function(x) if (length(x) >= 6) x[6] else NA),
    family  = sapply(split, function(x) if (length(x) >= 7) x[7] else NA),
    length = contig_lengths
  ) %>% 
  filter(!is.na(contig_ID))


plv_df <- read_sheet("https://docs.google.com/spreadsheets/d/113hSsqFV73bfdHTs5WoFFQU1DvnAcV5kPYmtanhROsY/edit?gid=0#gid=0", sheet = "PLVs")



# load tree data ----------------------------------------------------------
tree <- read.tree("intermediate/treebuilding/all_public_and_my_own_MCP_proteins_cleaned.aln.treefile")

tree_data <- data.table(tip_label = tree$tip.label,
                        short_id = "",
                        long_id = "",
                        class = "Virophaviricetes",
                        order = "",
                        family = "",
                        match_id = "",
                        origin = FALSE,
                        length = 0)


for(i in 1:nrow(tree_data)){
  # zero: if this is NP_048787.1 then we do something
  if(tree_data$tip_label[i] == "NP_048787.1"){
    # this is the root!
    tree_data$short_id[i] <- "PBCV-1"
    tree_data$long_id[i] <- "Paramecium bursaria Chlorella virus 1"
    tree_data$class[i] <- "Megaviricetes"
    tree_data$order[i] <- "Algavirales"
    tree_data$family[i] <- "Phycodnaviridae"
    tree_data$origin[i] <- "NCBI"
    tree_data$length[i] <- 0
    next
  }
  
  # first, if its one of ours
  if(str_detect(tree_data$tip_label[i], "\\_plv\\_|\\_vph\\_")){
    # is it a plv?
    if(str_detect(tree_data$tip_label[i], "\\_plv\\_")){
      # its a plv, fill in fields!
      plv_contig_name <- str_extract(tree_data$tip_label[i], "^.*\\_plv")
      tree_data$short_id[i] <- plv_df$public_ID[plv_df$contig == plv_contig_name]
      tree_data$long_id[i] <- tree_data$short_id[i]
      tree_data$class[i] <- "Polinton-like virus"
      tree_data$order[i] <- "Polinton-like virus"
      tree_data$family[i] <- "Unclassified"
      tree_data$origin[i] <- "this study"
      tree_data$length[i] <- plv_df$length[plv_df$public_ID == tree_data$short_id[i]]
      # skip to next one
      next
      
    }else{
      # then its a vph!
      vph_contig_name <- str_extract(tree_data$tip_label[i], "^.*\\_vph")
      tree_data$short_id[i] <- vph_df$public_ID[vph_df$contig_ID == str_remove(vph_contig_name, "\\_vph")]
      tree_data$long_id[i] <- tree_data$short_id[i]
      tree_data$class[i] <- "Virophaviricetes"
      tree_data$order[i] <- "Unclassified"
      tree_data$family[i] <- "Unclassified"
      tree_data$origin[i] <- "this study"
      tree_data$length[i] <- vph_df$length[vph_df$public_ID == str_remove(tree_data$short_id[i], "\\_vph")]
      
      # skip to next one
      next
    }
    
  }else{
    # then its a reference genome
    if(str_detect(tree_data$tip_label[i] , "assembled")){
      # then its a IMG genome
      tree_data$match_id[i] <- str_remove(str_remove(tree_data$tip_label[i], ".*assembled\\_"), "\\_\\d*$")
      tree_data$short_id[i] <- reference_df$short_id[str_detect(reference_df$long_id, tree_data$match_id[i])]
      tree_data$long_id[i] <- reference_df$long_id[str_detect(reference_df$long_id, tree_data$match_id[i])]
      tree_data$class[i] <- reference_df$class[match(tree_data$short_id[i], reference_df$short_id)]
      tree_data$order[i] <- reference_df$order[match(tree_data$short_id[i], reference_df$short_id)]
      tree_data$family[i] <- reference_df$family[match(tree_data$short_id[i], reference_df$short_id)]
      tree_data$origin[i] <- "IMG/VR3"
      tree_data$length[i] <- reference_df$length[match(tree_data$short_id[i], reference_df$short_id)]
      
      # skip to next one
      next
    }else{
      # then its either ncbi or hackl2021
      tree_data$match_id[i] <- str_remove(tree_data$tip_label[i], "\\_\\d*$")
      tree_data$long_id[i] <- tree_data$match_id[i]
      tree_data$short_id[i] <- reference_df$short_id[match(tree_data$match_id[i], reference_df$long_id)]
      tree_data$class[i] <- reference_df$class[match(tree_data$match_id[i], reference_df$long_id)]
      tree_data$order[i] <- reference_df$order[match(tree_data$match_id[i], reference_df$long_id)]
      tree_data$family[i] <- reference_df$family[match(tree_data$match_id[i], reference_df$long_id)]
      tree_data$length[i] <- reference_df$length[match(tree_data$short_id[i], reference_df$short_id)]
      
      
      if(str_detect(tree_data$short_id[i], "^EMALE.*")){
        tree_data$origin[i] <- "Hackl 2021"
      }else{
        tree_data$origin[i] <- "NCBI"
      }
      
      # skip to next one
      next
    }
  }
}

# make our lables bold
tree_data <- tree_data %>%
  mutate(label_bold = case_when(
    origin == "this study" ~ 2,
    TRUE ~ 1
  ))

# highlight the lavidavirales group to the left and the mividavirales cluster on the right
tree_data <- tree_data %>%
  mutate(highlight = case_when(
    # lavidavirales
    tip_label == "Lyne_tig00028020-10-249450_vph_26_27_28_concat" ~ "lavida",
    tip_label == "Damh_tig00046628-10-92530_vph_8" ~ "lavida",
    tip_label == "Lyne_tig00044463-10-176290_vph_16_15_concat" ~ "lavida",
    tip_label == "Lyne_tig00033829-10-240680_vph_12" ~ "lavida",
    tip_label == "Damh_tig00014446-10-220150_vph_23" ~ "lavida",
    tip_label == "Hade_tig00086668-10-71420_vph_7" ~ "lavida",
    # mividavirales
    tip_label == "Fred_tig00089364-10-137080_vph_6" ~ "mivida",
    tip_label == "Damh_tig00018526-10-141480_vph_14" ~ "mivida",
    tip_label == "Aved_tig00056523-10-134420_vph_14" ~ "mivida",
    tip_label == "MG807318.2_12" ~ "mivida", # zamilon
    tip_label == "JN603370.1_20" ~ "mivida", # sputnik 
    
    TRUE ~ NA
  ))





# plotting ----------------------------------------------------------------

# tree_midpoint <- phangorn::midpoint(tree = tree) # I used this only once, to show that the PLVs are the correct outgroup
# rerooted_tree <- root(tree, outgroup = "Ejby_tig00023995-10-186920_plv_14", resolve.root = TRUE)

rerooted_tree <- root(tree, outgroup = "NP_048787.1", resolve.root = TRUE)

# # test something
# rerooted_tree$edge.length[1] <- 0
# rerooted_tree$edge
# 
# r <- as.data.table(rerooted_tree$edge.length)
# t <- as.data.table(tree$edge.length)


a <- ggtree(rerooted_tree, layout = "fan", open.angle = 180, right = F, linewidth = 0.5) %<+% tree_data
a + 
  geom_tiplab(
    aes(label = ifelse(grepl("^IMG", short_id), NA, short_id), 
        fontface = label_bold), 
    align = TRUE, 
    linesize = 0.3, 
    size = 2, 
    offset = 1
  ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip_label, fill = origin),
    color = "black", offset = 0.05, size = 0.3, pwidth = 0.15
  ) +
  scale_fill_manual(
    values = c(
      `this study` = "steelblue",
      `IMG/VR3` = "grey90",
      NCBI = "grey60",
      `Hackl 2021` = "grey30"),
    guide = guide_legend(title = "Source", 
                         keywidth=0.5,
                         keyheight=0.5))  +
  
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = tip_label, fill = order),
    color = "black", offset = 0.05, size = 0.3, pwidth = 0.15
  ) +
  scale_fill_manual(
    values = c(
      Lavidavirales = "#ffdc02",
      Mividavirales = "brown",
      Priklausovirales = "#e99c05",
      Unclassified = "grey"),
    guide = guide_legend(title = "Order", 
                         keywidth=0.5,
                         keyheight=0.5)) +
  geom_fruit(
    geom=geom_bar,
    mapping=aes(y=tip_label, x=length),
    fill = "grey",
    # alpha = 0.5,
    pwidth=0.1, 
    orientation="y", 
    stat="identity",
    color="black",
    linewidth = 0.3,
    axis.params=list(
      axis       = "x",
      text.size  = 0,
      hjust      = 0,
      vjust      = 0.0,
      nbreak     = 2
    ),
    grid.params=list()
  ) +
  theme(legend.background=element_rect(fill=NA),
        legend.text=element_text(size=8) ,
        legend.position=c(0.92, 0.95),
        plot.margin = unit(c(0.3,1,-13,1), "cm")) +
  geom_nodepoint(
    mapping = aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), # Or > 0.95 for LPP
    color = "black",
    size = 1.0,
    shape = 19  # A solid circle
  ) +
  geom_rootedge(rootedge = 2)

ggsave(plot = last_plot(), file = "final/trees/plv_vph_references.pdf", height = 6, width = 10)
ggsave(plot = last_plot(), file = "final/trees/plv_vph_references.svg", height = 6, width = 10)
ggsave(plot = last_plot(), file = "final/trees/plv_vph_references.png", height = 6, width = 10)











