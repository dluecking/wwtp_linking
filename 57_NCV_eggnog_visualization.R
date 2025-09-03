# Author: dlu @ veelab
# Version: 2025-09-03

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(cowplot)
library(patchwork)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# helperfile data ---------------------------------------------------------

COG_df <- fread("helperfiles/cog_categories.tsv", sep = ",", header = F, )
names(COG_df) <- c("code", "description")


# load data ---------------------------------------------------------------

output_files <- list.files("intermediate/eggnog", pattern = "annotations", full.names = TRUE)
eggnog_df <- data.table()

for(file in output_files){
  # read file, skipe first 4 and last 3 rows
  tmp_df <- fread(file, fill = T, sep = "\t", skip = 4)
  tmp_df <- head(tmp_df, -3)
  
  eggnog_df <- rbind(eggnog_df, tmp_df)
}
# add description
eggnog_df$desc_short <- COG_df$description[match(eggnog_df$COG_category, COG_df$code)]



# data wrangling for plot -------------------------------------------------
overview_df <- data.table(type = c("function known", "function unknown"),
                          count = c(0, 0))

overview_df$count[overview_df$type == "function unknown"] <- eggnog_df %>% 
  filter(COG_category %in% c("R", "S", "RS", "SR", NA, "-")) %>% 
  nrow() %>% 
  unlist()

overview_df$count[overview_df$type == "function known"] <- nrow(eggnog_df) - eggnog_df %>% 
  filter(COG_category %in% c("R", "S", "RS", "SR", NA, "-")) %>% 
  nrow() %>% 
  unlist()



# visualization -----------------------------------------------------------

# first OVERVIEW plot
p1 <- ggplot(overview_df, aes(x = "", y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(
    values = c("function known" = "steelblue", "function unknown" = "grey"),
    labels = c("function known" = "Function known", "function unknown" = "Function unknown (R, S)")
  ) +
  labs(x = NULL, y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  theme_cowplot() +
  theme(axis.ticks.x = element_blank(),
        legend.position = "bottom")

 

# order the description
factor_levels <- c("Translation, ribosomal structure and biogenesis", "RNA processing and modification", "Transcription", "Replication, recombination and repair", "Chromatin structure and dynamics", "Cell cycle control, cell division, chromosome partitioning", "Nuclear structure", "Defense mechanisms", "Signal transduction mechanisms", "Cell wall/membrane/envelope biogenesis", "Cell motility", "Cytoskeleton", "Extracellular structures", "Intracellular trafficking, secretion, and vesicular transport", "Posttranslational modification, protein turnover, chaperones", "Energy production and conversion", "Carbohydrate transport and metabolism", "Amino acid transport and metabolism", "Nucleotide transport and metabolism", "Coenzyme transport and metabolism", "Lipid transport and metabolism", "Inorganic ion transport and metabolism", "Secondary metabolites biosynthesis, transport and catabolism")

# Convert the column to a factor with the specified levels
p2 <- ggplot(eggnog_df %>%
               mutate(desc_short = factor(desc_short, levels = factor_levels)) %>% 
               filter(!COG_category %in% c("R", "S", "RS", "SR", "-", "NA")) %>%
               filter(!is.na(desc_short)),
             aes(x = desc_short)) +
  geom_bar(color = "black", fill = "steelblue") +
  ylab("Count") +
  xlab(NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

p_combined <- p1 + p2 + plot_layout(widths = c(1, 6))  
p_combined

H <- 5
W <- 8
ggsave(plot = p_combined, file = "final/NCV_eggnog_barplot.png", height = H, width = W)
ggsave(plot = p_combined, file = "final/NCV_eggnog_barplot.pdf", height = H, width = W)
ggsave(plot = p_combined, file = "final/NCV_eggnog_barplot.svg", height = H, width = W)

