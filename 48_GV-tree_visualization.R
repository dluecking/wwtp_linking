# Author: dlu @ veelab
# Version: 2025-07-22

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

GV_info <- read_sheet("https://docs.google.com/spreadsheets/d/1QLNiqSt0XOS4xVPAeZAppwVjjjPIKdEE6w6f2_Qm55c/edit?gid=1228834474#gid=1228834474", 
                      sheet = "Final GVs overview") %>% 
  filter(!is.na(sample))

GV_info$tip_label <- GV_info$shortname


nuphylo_GV_data <- fread("helperfiles/NuPhylo_itol_colors.csv")


# add info to tree_data ---------------------------------------------------

# load tree first
aster_tree <- read.tree("intermediate/GV_tree/combined_nuphylo_astral.tre")

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
  left_join(GV_info %>% select(tip_label, personal_assessment_order, public_ID),
            by = "tip_label") %>%
  left_join(nuphylo_GV_data %>% select(Genome, Order),
            by = c("tip_label" = "Genome")) %>%
  mutate(
    tax_order = case_when(
      source == "this study" ~ personal_assessment_order,
      source == "reference" ~ Order
    )
  ) %>% 
  select(-Order, -personal_assessment_order)

# add the yaravirales one
tree_data$tax_order[tree_data$tip_label == "lcl|MT293574.1_prot_QKE44413.1"] <- "Yaravirales"

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
      tip_label == "Poxviridae_AF198100_Fowlpox_virus" ~ "Poxvirus Fowlpox (AF198100)",
      tip_label == "Phycodnaviridae_KF481685_Emiliania_huxleyi_virus_18" ~ "Emiliania_huxleyi_virus (KF481685)",
      tip_label == "Asfarviridae_KU702951_Faustovirus_strain_D6" ~ "Faustovirus D6 (KU702951)",
      tip_label == "Pithoviridae_MK072498_Solumvirus_sp_clone_Solumvirus_1" ~ "Pithovirus Solumnvirus 1(MK072498)",
      tip_label == "Iridoviridae_JQ724856_European_sheatfish_virus" ~ "Irodivirus Eur. Sheatfish V. (JQ724856)",
      tip_label == "Ascoviridae_KJ755191_Heliothis_virescens_ascovirus_3f_isolate_LD135790" ~ "Ascovirus Helithis virescens 3f",
      tip_label == "Marseilleviridae_HQ113105_Lausannevirus_isolate_7715" ~ "Lausannevirus 7715 (HQ113105)",
      tip_label == "Mimiviridae_ChoanoV1" ~ "Mimiviridae ChoanoV1",
      tip_label == "Mimiviridae_MF405918_Tupanvirus_deep_ocean" ~ "Tupanvirus deep ocean (MF405918)",
      tip_label == "Mimiviridae_KY684085_Indivirus_ILV1_Indivirus_1" ~ "Indivirus ILV1 (KY684085)",
      tip_label == "lcl|MT293574.1_prot_QKE44413.1" ~ "Yaravirus brasil. (MT293574.1)",
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

# first re-root to poxvirius (we dont have pokkes)
aster_tree <- ape::root(aster_tree, outgroup = "Poxviridae_AF198100_Fowlpox_virus", edgelabel = TRUE)

# midpoint rooting, only used once:
midpoint_tree <- phangorn::midpoint(aster_tree)

tree <- ggtree(aster_tree, layout = "fan", open.angle = 180)
#  flip(41, 186) # this was supposed to flip the one clade where the flipping would help, but it fucks up the tree...

# colors
paletteer_d("dichromat::BluetoOrangeRed_14")

tree %<+% tree_data + 
  geom_tiplab(mapping = aes(label = short_id, 
                            colour = show_in_tree, 
                            fontface = label_bold),
              linesize = 0.25,
              align = TRUE, 
              size = 1.55,
              hjust = 0,
              offset = 0.7) +
  # first order
  geom_fruit(geom = geom_tile,
             mapping = aes(y = tip_label, fill = tax_order),
             color = "black", offset = 0.02, size = 0.3) +
  scale_fill_manual(guide = guide_legend(title = "Order", ncol = 1),
                    values = c(
                      "Yaravirales"     = "#007A99FF",  # teal
                      "Asfuvirales"     = "#00AACCFF",  # orange
                      "Pimascovirales"  = "#66F0FFFF",  # purple
                      "Pandoravirales"  = "#CCFDFFFF",  # magenta
                      "Imitervirales"   = "#D9AF98FF",  # green
                      "Chitovirales"    = "#CC9B7AFF",  # mustard
                      "Proculvirales"   = "#F2DACEFF",  # brown
                      "Algavirales"     = "#662F00",  # blue
                      "unknown"         = "grey80"    # light grey
                    ),
                    na.translate = FALSE) +
  new_scale_fill() +
  # then completeness
  geom_fruit(geom = geom_tile,
             mapping = aes(y = tip_label, fill = completeness),
             color = "black", offset = 0.03, size = 0.3) +
  scale_fill_manual(
    values = c("complete" = "#333333", 
               "likely complete" = "#666666", 
               "likely incomplete" = "#999999", 
               "incomplete" = "#E5E5E5",
               "not assessed" = "white"),
    na.value = "white",
    na.translate = FALSE,  # suppress legend item and display for NA
    guide = guide_legend(title = "Completeness")
  ) +
  # Point showing circularity
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = tip_label, shape = circularity),
    fill = "black",
    alpha = 0.95,
    size = 1.2, offset = 0.03
  ) +
  # then labels
  scale_shape_manual(
    values = c("Y" = 1, "N" = 3),  # 21 = circle, 22 = square
    na.translate = FALSE,  # suppress legend item and display for NA,
    guide = guide_legend(title = "Circularity"),
    labels = c("Y" = "circular", "N" = "linear")
  ) +

  # scale_linetype_manual(values = c("blank" = "dotted", "aa" = "aa")) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "black"), guide = "none") +
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


ggsave(plot = last_plot(), file = "final/trees/GV_astral_w_references.pdf", height = 8, width = 8)
ggsave(plot = last_plot(), file = "final/trees/GV_astral_w_references.svg", height = 8, width = 8)



# # stats on different nuphylo trees ----------------------------------------
# 
# SINGLE_TREE_DIR <- list.dirs("intermediate/GV_tree/nuhpylo_output", recursive = FALSE)
# 
# trees <- lapply(paste0(SINGLE_TREE_DIR, "/allseqs.nwk"), read.tree)
# names(trees) <- basename(SINGLE_TREE_DIR)
# 
# library(TreeDist)
# 
# kc_distances <- KendallColijn(trees)






