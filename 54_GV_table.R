# Author: dlu @ veelab
# Version: 2025-08-26

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyverse)
library(gt)
library(gtExtras)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# # prep array stuff, do this once ------------------------------------------

# this gives us CRISPR spacer info
# c <- as.data.table(getName(read.fasta("intermediate/minced/minced_results_spacers_filtered_lc.fasta")))
# c$contig <- str_remove(c$V1, "\\_CRISPR.*$")
# 
# a <- as.data.table(table(c$contig)) 
# a <- a %>% 
#   filter(!str_ends(V1, "\\_lc$"))


# this is tRNA stuff
tRNA_df <- data.table()

for(file in list.files("intermediate/tRNA_prediction", full.names = TRUE)){
  if (length(readLines(file)) > 2) {
    tmp_df <- fread(file, skip = 2, header = F, sep = "\t")
  } else {
    next 
  }
  tmp_df$NCV <- str_remove(basename(file), "\\.trnas\\.txt$")
  tmp_df$tRNA <- str_remove(str_extract(tmp_df$V1, "tRNA-[A-z]*"), "tRNA-")
  tmp_df$start <- str_extract(tmp_df$V1, "(?<=\\[)\\d+")
  tmp_df$end <- str_extract(tmp_df$V1, "\\d+(?=\\])")
  tmp_df$anticodon <- str_extract(tmp_df$V3, "(?<=\\()[A-z]+(?=\\))")
  tmp_df$intein <- str_extract(tmp_df$V3, "i\\(\\d+,\\d+\\)")
  
  tmp_df <- tmp_df %>% 
    select(NCV, tRNA, start, end, anticodon, intein)
  tRNA_df <- rbind(tRNA_df, tmp_df)
}


# load crispr cas data ----------------------------------------------------

crispr_df <- fread("intermediate/CRISPR/cassette/HMM2019_cassettes.csv")
crispr_df$contig_id <- str_remove(crispr_df$V1, "\\_\\d+\\_ID.*$")


# load general data -------------------------------------------------------

# URL of the Google Sheet
sheet_url <- "https://docs.google.com/spreadsheets/d/1QLNiqSt0XOS4xVPAeZAppwVjjjPIKdEE6w6f2_Qm55c"

# Read the specified sheet and convert to data.table
gv_data <- read_sheet(sheet_url, sheet = "Final GVs overview")

gv_data <- gv_data %>% 
  filter(completeness %in% c("complete", "likely complete")) %>% 
  select(public_ID, shortname, sample, length, gc, personal_assessment_order, circular, completeness, ORFs, ncldv_hits, `tRNA (aragorn)`, padloc, crispr_array)

# fill in CRISPR data
gv_data$crispr_cas <- ""
for(i in 1:nrow(gv_data)){
  tmp_df <- crispr_df %>% 
    filter(contig_id == gv_data$shortname[i])
  
  if(nrow(tmp_df) >= 1){
    cas_genes <- tmp_df %>% 
      select(annotation) %>% 
      unlist() %>% 
      unique() %>% 
      paste(collapse = ", ")
  }else{
    cas_genes <- "-"
  }
  gv_data$crispr_cas[i] <- cas_genes
}



# fill in tRNA data to gv_data
gv_data$tRNA_list <- ""
gv_data$tRNA <- 0

for(i in 1:nrow(gv_data)){
  tmp_df <- tRNA_df %>% 
    filter(NCV == gv_data$shortname[i])
  
  if(nrow(tmp_df) == 0){
    next
  }
  
  gv_data$tRNA[i] <- nrow(tmp_df)
  gv_data$tRNA_list[i] <- paste(unique(tmp_df$tRNA[tmp_df$tRNA != ""]), collapse = ", ")
}


# create table ------------------------------------------------------------
# Corrected create table code with proper function order
gv_table <- gv_data %>%
  arrange(personal_assessment_order, completeness, public_ID) %>%
  mutate(length_plot = length) %>%
  select(
    public_ID, personal_assessment_order, completeness,
    length, length_plot, gc, ORFs, circular, ncldv_hits, padloc, 
    crispr_array, crispr_cas,
    tRNA, tRNA_list
  ) %>%
  gt() %>%
  
  # --- Column Visualizations ---
  gt_plt_bar(column = length_plot, color = "steelblue", keep_column = FALSE) %>%
  
  # move the bar column right next to Length (bp)
  cols_move(
    columns = length_plot,
    after = length
  ) %>%
  
  # move the crispr cas after the cas array
  cols_move(
    columns = crispr_cas,
    after = crispr_array
  ) %>%
  
  data_color(
    columns = gc,
    palette = c("lightblue", "steelblue", "darkblue"),
    domain = c(20, 60)
  ) %>%
  
  text_transform(
    locations = cells_body(columns = circular),
    fn = function(x) ifelse(x == "Y", "&#128902;", "&#9644;")
  ) %>%
  cols_align(align = "center", columns = circular) %>%
  
  gt_badge(
    column = completeness,
    palette = c("complete" = "darkgreen", "likely complete" = "lightgreen")
  ) %>%
  
  # --- General Formatting ---
  fmt_number(columns = length, decimals = 0) %>%
  
  tab_header(
    title = md("**Nucleocytoviruses identified in this study**"),
    subtitle = "Key genomic features of complete and likely complete NCVs."
  ) %>%
  
  cols_label(
    public_ID = "Genome ID",
    personal_assessment_order = "Predicted Order",
    completeness = "Completeness",
    length = "Length (bp)",
    length_plot = "Length (plot)",
    gc = "GC-content (%)",
    ORFs = "ORFs",
    circular = "Topology",
    ncldv_hits = "Marker Genes",
    tRNA = "# tRNAs",
    tRNA_list = "tRNA sequences",
    padloc = "Defensive System",
    crispr_array = "# CRISPR Spacers",
    crispr_cas = "CRISPR-Cas Genes"
  ) %>%
  
  fmt_missing(columns = everything(), missing_text = "-") %>%
  gt_theme_nytimes() %>%
  tab_style(
    style = cell_text(color = "#333333", weight = "bold"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(3)
  )


gv_table


gtsave(gv_table, filename = "final/ncv_table.html")
