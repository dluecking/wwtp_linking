# Author: dlu @ veelab
# Version: 2025-07-15

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(gggenomes)
library(googlesheets4)
library(gggenes)
library(paletteer)



# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load google sheet:
sheet_url <- "https://docs.google.com/spreadsheets/d/113hSsqFV73bfdHTs5WoFFQU1DvnAcV5kPYmtanhROsY/edit?usp=sharing"
vph_info <- read_sheet(sheet_url, sheet = "virophages") 

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

parse_prodigal_proteins_to_dataframe <- function(prodigal_file){
  tmp_fasta <- seqinr::read.fasta(prodigal_file)
  tmp_annotations <- unlist(seqinr::getAnnot(tmp_fasta))
  
  tmp_df <- data.table(
    prodigal_string = tmp_annotations,
    molecule = "",
    gene = "",
    start = 0,
    end = 0,
    strand = "",
    orientation = 0,
    gc_content = 0
  )
  
  tmp_df$molecule <- str_remove(str_remove(tmp_df$prodigal_string,"\\_\\d.*"), "^>")
  tmp_df$gene <- str_remove(str_remove(tmp_df$prodigal_string, " #.*"), "^>")
  tmp_df$start <- as.numeric(sapply(str_split(tmp_df$prodigal_string, " # "), `[`, 2))
  tmp_df$end <- as.numeric(sapply(str_split(tmp_df$prodigal_string, " # "), `[`, 3))
  tmp_df$strand <- as.numeric(sapply(str_split(tmp_df$prodigal_string, " # "), `[`, 4))
  tmp_df$gc_content <- as.numeric(sub(".*gc_cont=([0-9.]+)$", "\\1", sapply(str_split(tmp_df$prodigal_string, " # "), `[`, 5)))
  
  tmp_df <- tmp_df %>% 
    mutate(
      strand = case_when(
        strand == -1 ~ "reverse",
        TRUE ~ "forward"
      ),
      orientation = case_when(
        strand == "reverse" ~ 1,
        TRUE ~ 0
      )
    ) %>% 
    select(-prodigal_string)
  return(tmp_df)
}

# load data ---------------------------------------------------------------
# load hmm out if no local copy of df exists

if(!file.exists("local_data_storage/vph_hmm_out.tsv")){
  hmm_out <- data.table()
  
  for(file in list.files(path = "intermediate/hmm_out", pattern = "\\.tbl")){
    tmp_dt <- parseDomtblout(paste0("intermediate/hmm_out/", file))
    tmp_dt$sample <- str_remove(file, "\\_.*$")
    
    hmm_out <- rbind(hmm_out, tmp_dt)
  }
  rm(tmp_dt)
  
  hmm_out <- hmm_out %>% 
    filter(!is.na(query)) %>% 
    filter(score >= 50) %>% 
    filter(query != "ALL_PLV")
  
  
  hmm_out$contig <- str_remove(hmm_out$target, "\\_\\d*$")
  hmm_out$query_type <- str_remove(hmm_out$query, "\\_\\d*$")
  hmm_out$gene_id <- paste0(hmm_out$sample, "_", hmm_out$target)
  hmm_out$gene_id <- str_replace(hmm_out$gene_id, "_(?!.*_)", "_vph_")
  
  fwrite(hmm_out, file = "local_data_storage/vph_hmm_out.tsv")
}else{
  hmm_out <- fread("local_data_storage/vph_hmm_out.tsv")
}

# load prodigal orfs if no local copy exists
if(!file.exists("local_data_storage/vph_gene_df.tsv")){
  prodigal_files <- list.files("intermediate/proteins/vph", full.names = TRUE)
  
  gene_df <- data.table()
  
  for(file in prodigal_files){
    tmp_df <- parse_prodigal_proteins_to_dataframe(file)
    gene_df <- rbind(gene_df, tmp_df)
  }
  rm(tmp_df)
  
  # combine gene info with the hmm_out
  gene_df$annotation <- hmm_out$query[match(gene_df$gene, hmm_out$gene_id)]
  gene_df$annotation[is.na(gene_df$annotation)] <- "hypothetical"
  
  fwrite(gene_df, "local_data_storage/vph_gene_df.tsv")
  
}else{
  gene_df <- fread("local_data_storage/vph_gene_df.tsv")
}

# add introscan annotation, load local if exists
if(!file.exists("local_data_storage/vph_interpro_out.tsv")){
  
  interpro_files <- list.files("intermediate/annotations/interpro/vph", pattern = ".tsv$", full.names = TRUE)
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
  
  fwrite(interpro_out, "local_data_storage/vph_interpro_out.tsv")
}else{
  interpro_out <- fread("local_data_storage/vph_interpro_out.tsv")
}

# add interpro annotation but only if not already known
for(i in 1:nrow(gene_df)){
  # if not empty interpro annotation AND the currecnt annotation is "hypothetical"
  if(length(interpro_out$V6[interpro_out$V1 == gene_df$gene[i]]) > 0 & gene_df$annotation[i] == "hypothetical"){
    gene_df$annotation[i] <- interpro_out$V6[interpro_out$V1 == gene_df$gene[i]]
  }
}

# simplify annotations
gene_df <- gene_df %>% 
  mutate(annotation_short = case_when(
    annotation == "DNA/RNA polymerases" ~ "DNA Pol",
    annotation == "Palm domain of DNA polymerase" ~ "DNA Pol",
    annotation == "His-Me finger endonucleases" ~ "His-Me ENase",
    annotation == "Major capsid protein V20 C-terminal domain" ~ "MCP C-terminal domain",
    annotation == "P-loop containing nucleoside triphosphate hydrolases" ~ "P-loop NTPase",
    annotation == "S-adenosyl-L-methionine-dependent methyltransferases" ~ "SAM-dependent_MTases_sf",
    annotation == "Ribonuclease H-like" ~ "RNaseH-like sf",
    annotation == "Sputnik minor capsid protein V18/19" ~ "mCP V18/19",
    annotation == "Sputnik virophage major capsid protein 1st domain" ~ "MCP 1st domain",
    TRUE ~ annotation
    )
  )


# adjust orientation based on Penton
gene_df_reversed <- gene_df

for(seq in unique(gene_df_reversed$molecule)){
  PENTON_STRAND <- gene_df_reversed %>% 
    filter(molecule == seq) %>% 
    filter(annotation_short == "Penton_1") %>% 
    sample_n(1) %>% 
    pull(strand)
  
  # only if MCP is reverse
  if(PENTON_STRAND == "reverse"){
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


# visualization of vphs ---------------------------------------------------
# all annotations which are only in a single vph get "other"

gene_df_reversed <- gene_df_reversed %>%
  group_by(annotation_short) %>%
  mutate(count = n()) %>% # Calculate the count for each annotation_short
  ungroup() %>%
  mutate(annotation_short = ifelse(count < 2, "other", annotation_short)) %>%
  select(-count)

# 1. Add the "_vph" suffix to the contig_ID column in vph_info
vph_info <- vph_info %>%
  mutate(contig_ID = paste0(contig_ID, "_vph"))

# Now, create the mapping from the corrected data frame
id_map <- setNames(vph_info$public_ID, vph_info$contig_ID)

# 2. Add the new public IDs to the main data frame
gene_df_reversed <- gene_df_reversed %>%
  mutate(
    # Create a base name by removing the " (rev)" suffix for matching
    base_molecule = str_remove(molecule, " \\(rev\\)"),
    # Map the base name to the new public ID
    plot_id_base = id_map[base_molecule],
    # Conditionally add " (rev)" back if the original molecule was reversed
    plot_id = ifelse(str_detect(molecule, " \\(rev\\)"),
                     paste0(plot_id_base, " (rev)"),
                     plot_id_base)
  )

# 3. Define the original order (as before)
desired_molecule_order <- c(
  "Hade_tig00086668-10-71420_vph (rev)",
  "Damh_tig00014446-10-220150_vph", # cluster 2
  "Lyne_tig00033829-10-240680_vph (rev)",
  "Lyne_tig00044463-10-176290_vph",
  "Damh_tig00046628-10-92530_vph (rev)",
  "Lyne_tig00028020-10-249450_vph",
  "Vibo_tig00073545-10-42870_vph",
  "AalE_tig00021708-10-192480_vph",
  "Fred_tig00089364-10-137080_vph (rev)",
  "Damh_tig00018526-10-141480_vph", # cluster 3
  "Aved_tig00056523-10-134420_vph",
  

  "Fred_tig00051270-10-161640_vph (rev)",
  "Skiv_tig00138093-10-67990_vph (rev)",
  "Viby_tig00043866-10-124870_vph (rev)",
  "Aved_tig00303955-10-54120_vph" # cluster 1
)

# 4. Create the corresponding order for the new public IDs
base_order <- str_remove(desired_molecule_order, " \\(rev\\)")
desired_plot_id_order_base <- id_map[base_order]
desired_plot_id_order <- ifelse(str_detect(desired_molecule_order, " \\(rev\\)"),
                                paste0(desired_plot_id_order_base, " (rev)"),
                                desired_plot_id_order_base)

# 5. Apply the new ordering to the plot_id column by making it a factor
gene_df_reversed$plot_id <- factor(gene_df_reversed$plot_id, levels = desired_plot_id_order)


# 6. Create dummies for nice alignment using the new plot_id
dummies <- make_alignment_dummies(
  gene_df_reversed %>% select(plot_id, gene, start, end, annotation_short),
  aes(xmin = start, xmax = end, y = plot_id, id = annotation_short),
  on = "Penton_1"
)


# visualization of vphs ---------------------------------------------------
ggplot(gene_df_reversed, aes(xmin = start, xmax = end, y = plot_id, fill = annotation_short)) +
  geom_gene_arrow(aes(forward = orientation), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), size = 0.5) +
  geom_blank(data = dummies) +
  facet_wrap(~ plot_id, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(
    # vph only genes (core genes)
    MCP_1 = "goldenrod1",
    ATPase_1 = "#529985",
    Penton_1 = "#689D72",
    PRO_1 = "#9EB258",
    
    # then vph only genes 
    `Integrase core domain` = "#C26B51", 
    `SET domain` = "#D78D50",
    `SGNH hydrolase` = "#E5A94E",
    `alpha/beta-Hydrolases` = "#E7C64A",
    `Concanavalin A-like lectins/glucanases` = "#DEC895",
    other = "white",                   # Light gray for miscellaneous
    
    # shared:
    `DNA Pol` = "#333333",
    `His-Me ENase` = "#666666",
    `P-loop NTPase` = "#AAAAAA",
    `RNaseH-like sf` = "#DDDDDD",
    hypothetical = "white"
    

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
    legend.position = "none",
    axis.text.y = element_text(face = "bold")
  )

H = 4
W = 7
ggsave(plot = last_plot(), file = "final/vph_genome_map.png", height = H, width = W)
ggsave(plot = last_plot(), file = "final/vph_genome_map.svg", height = H, width = W)
ggsave(plot = last_plot(), file = "final/vph_genome_map.pdf", height = H, width = W)



palette_used <- paletteer_d("ggthemes::Tableau_20")

