# Author: dlu @ veelab
# Version: 2025-07-11

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(gggenomes)
library(googlesheets4)
library(gggenes)



# set wogggenomes# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# functions ---------------------------------------------------------------

parseDomtblout <- function(file_path){
  con <- file(file_path, "r")
  lines <- c("")
  while (TRUE) {
    l <- readLines(con, n = 1)
    
    # Break the loop if the end of the file is reached
    if (length(l) == 0) {
      break
    }
    # append non-comment lines to lines
    if(!startsWith(l, "#")){
      lines <- c(lines, l)
    }
  }
  # Close the file connection
  close(con)
  
  # now take the lines, and manipulate each of them, so we return a nice dt
  lines <- as.data.table(lines) %>% 
    filter(lines != "")
  
  dt <- data.table(target = as.character(),
                   t_acc = as.character(),
                   query = as.character(),
                   q_acc = as.character(),
                   evalue = as.numeric(),
                   score = as.numeric())
  for(i in 1:nrow(lines)){
    tmp <- unlist(str_split(lines[i], " "))
    tmp <- tmp[tmp!=""]
    # select first 6 elements
    tmp <- tmp[1:6]
    tmp_dt <- data.table(target = tmp[1],
                         t_acc = tmp[2],
                         query = tmp[3],
                         q_acc = tmp[4],
                         evalue = as.numeric(tmp[5]),
                         score = as.numeric(tmp[6]))
    dt <- rbind(dt, tmp_dt)
  }
  return(dt)
}

parse_prokka_output <- function(file_path) {
  raw_output <- readLines(file_path) %>%
    str_trim() %>%
    str_subset("^$", negate = TRUE)
  
  parsed_features <- list()
  current_seq_id <- NA
  
  for (line in raw_output) {
    if (str_detect(line, "^Feature ")) {
      current_seq_id <- str_remove(line, "^Feature ")
    } else if (str_detect(line, "^\\d+")) {
      parts <- str_split(line, "\\s+", n = 3, simplify = TRUE) %>% str_trim()
      start_raw <- as.numeric(parts[1])
      end_raw <- as.numeric(parts[2])
      type <- parts[3]
      
      strand <- if (start_raw < end_raw) "+" else "-"
      start <- min(start_raw, end_raw)
      end <- max(start_raw, end_raw)
      
      parsed_features[[length(parsed_features) + 1]] <- list(
        file_id = basename(file_path),
        seq_id = current_seq_id,
        start = start,
        end = end,
        strand = strand,
        type = type,
        locus_tag = NA_character_,
        inference = NA_character_,
        product = NA_character_,
        introns = list(NULL),
        parent_ids = list(NULL),
        source = "Prokka",
        score = NA_real_,
        phase = NA_integer_,
        width = end - start + 1,
        gc_content = NA_real_,
        name = NA_character_,
        Note = NA_character_
      )
    } else if (str_detect(line, "^\\s+")) {
      last_idx <- length(parsed_features)
      if (last_idx == 0) next
      
      parts <- str_split(str_trim(line), "\\s+", n = 2, simplify = TRUE)
      key <- tolower(parts[1])
      value <- parts[2]
      
      if (key == "locus_tag") {
        parsed_features[[last_idx]]$locus_tag <- value
      } else if (key == "inference") {
        parsed_features[[last_idx]]$inference <- value
      } else if (key == "product") {
        parsed_features[[last_idx]]$product <- value
      }
    }
  }
  
  df <- bind_rows(parsed_features) %>%
    # Ensure seq_id is never NA, providing a fallback if the 'Feature' line was missing.
    mutate(
      seq_id = coalesce(seq_id, str_remove(basename(file_id), "\\.tbl$"))
    ) %>%
    # Generate sequential IDs per contig and per feature type for robust fallback.
    group_by(seq_id, type) %>%
    mutate(
      sequential_id_for_type = row_number()
    ) %>%
    ungroup() %>%
    mutate(
      # Prioritize locus_tag. If missing, use seq_id_type_sequential_number.
      feat_id = coalesce(locus_tag, paste0(seq_id, "_", sequential_id_for_type)),
      geom_id = feat_id, # geom_id typically matches feat_id
      # Name should also prioritize locus_tag, then product, then a generic type_sequential_id.
      name = coalesce(locus_tag, product, paste0(type, "_", sequential_id_for_type)),
      # Note can prioritize product or other fields as needed.
      Note = coalesce(product, Note)
    ) %>%
    select(
      file_id, seq_id, start, end, strand, type, feat_id, introns, parent_ids,
      source, score, phase, width, gc_content, name, Note, geom_id
    )
  
  return(df)
}


# load data ---------------------------------------------------------------

hmm_out <- data.table()

for(file in list.files(path = "intermediate/hmm_out", pattern = "\\.tbl")){
  tmp_dt <- parseDomtblout(paste0("intermediate/hmm_out/", file))
  tmp_dt$sample <- str_remove(file, "\\_.*$")
  
  hmm_out <- rbind(hmm_out, tmp_dt)
  
  
}
rm(tmp_dt)

hmm_out <- hmm_out %>% 
  filter(!is.na(query)) %>% 
  filter(score >= 50)


hmm_out$contig <- str_remove(hmm_out$target, "\\_\\d*$")
hmm_out$query_type <- str_remove(hmm_out$query, "\\_\\d*$")



# first PLV ---------------------------------------------------------------

sheet_url <- "https://docs.google.com/spreadsheets/d/113hSsqFV73bfdHTs5WoFFQU1DvnAcV5kPYmtanhROsY/edit?usp=sharing"
plv_info <- read_sheet(sheet_url, sheet = "PLVs") 
plv_hmm <- hmm_out %>% filter(query == "ALL_PLV")
plv_hmm$gene_id <- paste0(plv_hmm$sample, "_", plv_hmm$target)
plv_hmm$gene_id <- str_replace(plv_hmm$gene_id, "_(?!.*_)", "_plv_")

prokka_dirs <-list.dirs("intermediate/annotations/prokka/plv", recursive = F)
prokka_files <- paste0(prokka_dirs,"/", basename(prokka_dirs), ".tbl")

gene_df <- data.table()
for(file in prokka_files){
  tmp_df <- parse_prokka_output(file) %>% filter(type == "gene")
  gene_df <- rbind(gene_df, tmp_df)
}

# combine gene info with the hmm_out
gene_df$Note <- plv_hmm$query[match(gene_df$feat_id, plv_hmm$gene_id)]
gene_df$Note[gene_df$Note == "ALL_PLV"] <- "MCP"

# get this into gggenes format
gene_df <- gene_df %>% 
  mutate(molecule = seq_id,
         gene = feat_id,
         start = start,
         end = end,
         strand = strand,
         orientation = strand,
         annotation = Note
  ) %>% 
  select(molecule, gene, start, end, strand, orientation, annotation)

# fix orientation
gene_df <- gene_df %>% 
  mutate(
    strand = case_when(
      strand == "-" ~ "reverse",
      TRUE ~ "forward"),
    orientation = case_when(
      orientation == "-" ~ 1,
      TRUE ~ 0
    ),
    annotation = case_when(
      is.na(annotation) ~ "hypothetical",
      TRUE ~ annotation
    )
  )


# add introscan annotation
interpro_files <- list.files("intermediate/annotations/interpro/plv", pattern = ".tsv$", full.names = TRUE)
interpro_out <- rbindlist(lapply(interpro_files, fread))

# filter out unnecessary hits
interpro_out <- interpro_out %>% 
  filter(V9 != "-") %>% 
  mutate(evalue = as.numeric(V9)) %>% 
  filter(evalue <= 10^-5) %>% 
  filter(V6 != "-" & V13 != "-") %>% 
  filter(!str_detect(V6, "unknown function"))

# only get best hit per id based on evalue
interpro_out <- interpro_out %>% 
  group_by(V1) %>% 
  slice_min(V9, n = 1)

# add interpro annotation but only if not MCP
for(i in 1:nrow(gene_df)){
  # if not empty interpro annotation AND not MCP annotation already present
  if(length(interpro_out$V6[interpro_out$V1 == gene_df$gene[i]]) > 0 & gene_df$annotation[i] != "MCP"){
    gene_df$annotation[i] <- interpro_out$V6[interpro_out$V1 == gene_df$gene[i]]
  }
}

# simplyfy annotations

gene_df <- gene_df %>% 
  mutate(annotation_short = case_when(
    annotation == "hypothetical" ~ "hypothetical",
    annotation == "MCP" ~ "MCP",
    annotation == "DNA/RNA polymerases" ~ "DNA Pol",
    annotation == "DNA polymerase type B, organellar and viral" ~ "DNA Pol",
    annotation == "Palm domain of DNA polymerase" ~ "DNA Pol",
    annotation == "His-Me finger endonucleases" ~ "His-Me ENase",
    annotation == "GIY-YIG endonuclease" ~ "GIY-YIG ENase",
    annotation == "P-loop containing nucleoside triphosphate hydrolases" ~ "P-loop NTPase",
    annotation == "Poxvirus A32 protein" ~ "A32",
    annotation == "Ribonuclease H-like" ~ "RNaseH-like sf")
  )


# visualize PLV -----------------------------------------------------------

# Create a new data frame for the reversed coordinates
gene_df_reversed <- gene_df

for(seq in unique(gene_df_reversed$molecule)){
  MCP_STRAND <- gene_df_reversed %>% 
    filter(molecule == seq) %>% 
    filter(annotation_short == "MCP") %>% 
    pull(strand)
  
  # only if MCP is reverse
  if(MCP_STRAND == "reverse"){
    # total length of current molecule
    seq_length <- max(gene_df_reversed$end[gene_df_reversed$molecule == seq])
    
    # which indices are this for the gene_df?
    row_indices <- which(gene_df_reversed$molecule == seq)
    
    # store original coords
    old_start <- gene_df_reversed$start[row_indices]
    old_end <- gene_df_reversed$end[row_indices]
    
    # calc new coords
    gene_df_reversed$start[row_indices] <- seq_length - old_end
    gene_df_reversed$end[row_indices] <- seq_length - old_start
    
    # reverse strand
    old_strand <- gene_df_reversed$strand[row_indices]
    gene_df_reversed$strand[row_indices] <- ifelse(old_strand == "forward", "reverse", "forward")
    
    # flip orientation
    old_orientation <- gene_df_reversed$orientation[row_indices]
    gene_df_reversed$orientation[row_indices] <- ifelse(old_orientation == 0, 1, 0)
    
    # adjust the name accordingly
    gene_df_reversed$molecule[row_indices] <- paste0(gene_df_reversed$molecule[row_indices], " (rev)")
    
  }
}

# order needs to be set
# this is sorted by subclusters on top (1, 2, 3, then singles)
desired_molecule_order <- c(
  "AalW_tig00083928-10-131450_plv (rev)",
  "Rand_tig00054083-10-84480_plv (rev)",
  "Vibo_tig00019442-10-136670_plv (rev)",
  "Bjer_tig00027726-10-154320_plv (rev)",
  "Mari_tig00039850-10-191400_plv",
  "Ribe_tig00030249-10-170680_plv",
  "Aved_tig00048883-10-193450_plv (rev)",
  "Aved_tig00084897-10-150130_plv (rev)",
  "Ejby_tig00023995-10-186920_plv (rev)",
  "Lyne_tig00046032-10-146580_plv (rev)",
  "Lyne_tig00060056-10-166770_plv",
  "Rand_tig00055912-10-108660_plv",
  "Rand_tig00813462-10-110980_plv",
  "Vibo_tig00024931-10-107240_plv"
)
gene_df_reversed$molecule <- factor(gene_df_reversed$molecule, levels = desired_molecule_order)


# create dummies for nice alignment
dummies <- make_alignment_dummies(
  gene_df_reversed,
  aes(xmin = start, xmax = end, y = molecule, id = annotation_short),
  on = "MCP"
)

ggplot(gene_df_reversed, aes(xmin = start, xmax = end, y = molecule, fill = annotation_short)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(
    MCP = "hotpink",            # Stays hotpink - distinct and vibrant
    `DNA Pol` = "#014263",      # Deep, bold red - excellent for highlighted/important
    A32 = "#00ffba",            # Soft, light blue - subtle, not highlighted
    `GIY-YIG ENase` = "#ffa05f", # Lighter, softer purple
    `His-Me ENase` = "#fac55b",  # Deeper, richer purple - clearly related, but distinct
    `P-loop NTPase` = "#2dd2c0", # Warm, rich burnt orange
    `RNaseH-like sf` = "#fc8484", # Deeper, more brownish-orange - clearly related, but distinct
    hypothetical = "white"    # Very light grey, almost white - for background/hypothetical
  )) +
  theme_genes() +
  labs(y = "") +
  theme(
    axis.text.x = element_blank(),       # Removes the axis numbers/labels
    axis.ticks.x = element_blank(),      # Removes the tick marks
    axis.title.x = element_blank(),      # Removes the x-axis title (if any)
    panel.grid.major.x = element_blank(), # Removes major vertical grid lines
    panel.grid.minor.x = element_blank(),  # Removes minor vertical grid lines
    axis.line.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

ggsave(plot = last_plot(), file = "final/plv_genome_map.png", height = 4.5, width = 8)

