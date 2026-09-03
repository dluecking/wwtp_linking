# Author: dlu @ veelab
# Version: 2026-07-17
# Same visualization as REVIEW_03, but on the taxon-sampling astral tree
# (REVIEW_06/07/08 output) with yaraviruses removed - so no yara handling
# needed anymore.

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(ggtree)
library(ggtreeExtra)
library(googlesheets4)
library(ggnewscale)
library(paletteer)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load auxiliary data -----------------------------------------------------

GV_info <- read_sheet("https://docs.google.com/spreadsheets/d/159KaIRGjUGnN8uIJVsG42HIvWQ8MG2dysxNbNqw0VBg/edit?gid=697074736#gid=697074736",
                      sheet = "Table S1")

GV_info$tip_label <- GV_info$shortname


nuphylo_GV_data <- fread("helperfiles/NuPhylo_itol_colors.csv")
nuphylo_GV_data$Order[nuphylo_GV_data$Order == "Pandoravirales"] <- '"Pandoravirales"'

# add info to tree_data ---------------------------------------------------

# load tree first (taxon-sampling astral tree, no yaraviruses)
aster_tree <- read.tree("intermediate/GV_tree/REVIEW_taxon_sampling/astral/combined_taxon_sampling_astral.tre")

# create tree_data
tree_data <- data.table(tip_label = unique(c(GV_info$tip_label, aster_tree$tip.label)))


# make them bold if our study
tree_data <- tree_data %>%
  mutate(source = case_when(
    tip_label %in% GV_info$tip_label ~ "this study",
    TRUE ~ "reference"
  )) %>%
  mutate(label_bold = case_when(
    tip_label %in% GV_info$tip_label ~ 2,
    TRUE ~ 1
  ))



# but now lets add family information so we can color!
tree_data <- tree_data %>%
  left_join(GV_info %>% select(tip_label, order, public_ID),
            by = "tip_label") %>%
  left_join(nuphylo_GV_data %>% select(Genome, Order),
            by = c("tip_label" = "Genome")) %>%
  mutate(
    tax_order = case_when(
      source == "this study" ~ order,
      source == "reference" ~ Order
    )
  ) %>%
  select(-Order, -order)

# add info of our genomes
# circularity
tree_data$circularity <- GV_info$circular[match(tree_data$tip_label, GV_info$tip_label)]

# completeness estimate
tree_data$completeness <- GV_info$completeness[match(tree_data$tip_label, GV_info$tip_label)]
tree_data$completeness[is.na(tree_data$completeness)] <- "not assessed"


# short_id needs to be added: you are "" unless you are in my dataset, OR I manually want to higlight with a specific short_id
tree_data <- tree_data %>%
  mutate(
    short_id = case_when(
      tip_label %in% GV_info$tip_label ~ public_ID,
      tip_label == "Poxviridae_AF198100_Fowlpox_virus" ~ "ROOT: Poxvirus Fowlpox (AF198100)",
      tip_label == "Phycodnaviridae_KF481685_Emiliania_huxleyi_virus_18" ~ "Emiliania_huxleyi_virus (KF481685)",
      tip_label == "Asfarviridae_KJ614390_Faustovirus_strain_E12" ~ "Faustovirus E12 (KJ614390)",
      tip_label == "TARA_ARC_NCLDV_00085" ~ "TARA_ARC_NCLDV_00085",
      tip_label == "Iridoviridae_JQ724856_European_sheatfish_virus" ~ "Irodivirus Eur. Sheatfish V. (JQ724856)",
      tip_label == "Ascoviridae_KJ755191_Heliothis_virescens_ascovirus_3f_isolate_LD135790" ~ "Ascovirus Helithis virescens 3f",
      tip_label == "Marseilleviridae_HQ113105_Lausannevirus_isolate_7715" ~ "Lausannevirus 7715 (HQ113105)",
      tip_label == "Mimiviridae_ChoanoV1" ~ "Mimiviridae ChoanoV1",
      tip_label == "Mimiviridae_MF405918_Tupanvirus_deep_ocean" ~ "Tupanvirus deep ocean (MF405918)",
      tip_label == "Mimiviridae_KY684085_Indivirus_ILV1_Indivirus_1" ~ "Indivirus ILV1 (KY684085)",
      tip_label == "Phycodnaviridae_JX997184_Paramecium_bursaria_Chlorella_virus_OR0704_2_2" ~ "Paramecium bursaria Chlorella virus (OR0704)",
      tip_label == "Poxviridae_AF170722_Rabbit_fibroma_virus" ~ "Rabbit fibroma virus (AF170722)",
      TRUE ~ ""
    )
  )

# which labels to show in the tree?
tree_data <- tree_data %>%
  mutate(
    show_in_tree = case_when(
      short_id != "" ~ "TRUE",
      TRUE ~ "FALSE"
    )
  )

# which which line type do we want to have
tree_data <- tree_data %>%
  mutate(
    line_type = ifelse(short_id == "", "dotted", "aa")
  )
# which which line type do we want to have
tree_data <- tree_data %>%
  mutate(
    numeric_linesize = ifelse(short_id == "", "0", "0.05")
  )


# vis tree ----------------------------------------------------------------

# change procul to "procul"
tree_data$tax_order[tree_data$tax_order == "Proculvirales"] <- '"Proculvirales"'

# first re-root to poxvirius (we dont have pokkes) - same outgroup as REVIEW_03
aster_tree <- ape::root(aster_tree, outgroup = "Poxviridae_AF198100_Fowlpox_virus", edgelabel = TRUE)

tree <- ggtree(aster_tree, layout = "fan", open.angle = 90)

# colors
paletteer_d("dichromat::BluetoOrangeRed_14")


tree %<+% tree_data +
  geom_tiplab(mapping = aes(label = short_id,
                            colour = show_in_tree,
                            fontface = label_bold),
              linesize = 0.25,
              align = TRUE,
              size = 2.5,
              hjust = 0,
              offset = 1.5) +
  # first order
  geom_fruit(geom = geom_tile,
             mapping = aes(y = tip_label, fill = tax_order),
             color = "black", offset = 0.03, size = 0.3, pwidth = 0.4) +
  scale_fill_manual(guide = guide_legend(title = "Order", ncol = 1),
                    values = c(
                      "Asfuvirales"     = "#00AACCFF",  # orange
                      "Pimascovirales"  = "#00d6ed",  # purple
                      '"Pandoravirales"'  = "#CCFDFFFF",  # magenta
                      "Pandoravirales"  = "#CCFDFFFF",  # magenta
                      "Imitervirales"   = "#b4934b",  # green
                      "Chitovirales"    = "#CC9B7AFF",  # mustard
                      '"Proculvirales"'   = "#F2DACEFF",  # brown
                      "Algavirales"     = "#662F00",  # blue
                      "unknown"         = "grey80"    # light grey
                    ),
                    na.translate = FALSE) +
  new_scale_fill() +
  # then completeness
  geom_fruit(geom = geom_tile,
             mapping = aes(y = tip_label, fill = completeness),
             color = "black", offset = 0.03, size = 0.3, pwidth = 0.4) +
  scale_fill_manual(
    values = c("complete" = "grey20",      # Blue (high quality)
               "draft" = "grey80",         # Light blue (lower quality)
               "not assessed" = "white"),   # White (no assessment)
    na.value = "white",
    guide = guide_legend(title = "Completeness")
  ) +
  # Point showing circularity
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = tip_label, shape = circularity),
    color = "black",
    alpha = 0.95,
    size = 2.5, offset = 0.04, pwidth = 0.4
  ) +
  # then labels
  scale_shape_manual(
    values = c("Y" = 16, "N" = 4),  # 21 = circle, 22 = square
    na.translate = FALSE,  # suppress legend item and display for NA,
    guide = guide_legend(title = "Circularity"),
    labels = c("Y" = "circular", "N" = "linear")
  ) +
  # scale_linetype_manual(values = c("blank" = "dotted", "aa" = "aa")) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"), guide = "none") +
  theme_tree() +
  theme(legend.background=element_rect(fill=NA),
        legend.text=element_text(size=6),
        legend.title = element_text(size = 7, face = "bold"),
        legend.position=c(1, 0.28),
        legend.key.size = unit(0.25, "cm"),
        plot.margin = unit(c(2,2,0.5,0.8), "cm"),
        legend.spacing.y = unit(0.1, "cm")) +
  xlim_tree(c(NA, 5)) +
  geom_nodepoint(
    mapping = aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 0.75), # Or > 0.95 for LPP
    color = "black",
    size = 1.0,
    shape = 19
  )


ggsave(plot = last_plot(), file = "final/trees/REVIEW_taxon_sampling_GV_astral_w_references.pdf", height = 8, width = 8)
ggsave(plot = last_plot(), file = "final/trees/REVIEW_taxon_sampling_GV_astral_w_references.svg", height = 8, width = 8)
