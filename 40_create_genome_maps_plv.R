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

# but only if a local df verison of this one does not exist already...
if(!file.exists("local_data_storage/plv_hmm_out.tsv")){
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
  
  fwrite(hmm_out, file = "local_data_storage/plv_hmm_out.tsv")
}else{
  hmm_out <- fread("local_data_storage/plv_hmm_out.tsv")
}





# first PLV ---------------------------------------------------------------

sheet_url <- "https://docs.google.com/spreadsheets/d/113hSsqFV73bfdHTs5WoFFQU1DvnAcV5kPYmtanhROsY/edit?usp=sharing"
plv_info <- read_sheet(sheet_url, sheet = "PLVs") 

plv_hmm <- hmm_out %>% filter(query == "ALL_PLV")
plv_hmm$gene_id <- paste0(plv_hmm$sample, "_", plv_hmm$target)
plv_hmm$gene_id <- str_replace(plv_hmm$gene_id, "_(?!.*_)", "_plv_")


# read gene_df if no local copy exists
if(!file.exists("local_data_storage/plv_gene_df.tsv")){
  prodigal_files <- list.files("intermediate/proteins/plv", full.names = TRUE)
  
  gene_df <- data.table()
  
  for(file in prodigal_files){
    tmp_df <- parse_prodigal_proteins_to_dataframe(file)
    gene_df <- rbind(gene_df, tmp_df)
  }
  rm(tmp_df)
  
  # combine gene info with the hmm_out
  gene_df$annotation <- plv_hmm$query[match(gene_df$gene, plv_hmm$gene_id)]
  gene_df$annotation[gene_df$annotation == "ALL_PLV"] <- "MCP"
  gene_df$annotation[is.na(gene_df$annotation)] <- "hypothetical"
  
  fwrite(gene_df, "local_data_storage/plv_gene_df.tsv")
}else{
  gene_df <- fread("local_data_storage/plv_gene_df.tsv")
}


# add introscan annotation, read the interpro out, unless local copy exsists
if(!file.exists("local_data_storage/plv_interpro_out.tsv")){
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
  
  fwrite(interpro_out, "local_data_storage/plv_interpro_out.tsv")
}else{
  interpro_out <- fread("local_data_storage/plv_interpro_out.tsv")
}


# add interpro annotation but only if not MCP
for(i in 1:nrow(gene_df)){
  # if not empty interpro annotation AND not MCP annotation already present
  if(length(interpro_out$V6[interpro_out$V1 == gene_df$gene[i]]) > 0 & gene_df$annotation[i] != "MCP"){
    gene_df$annotation[i] <- interpro_out$V6[interpro_out$V1 == gene_df$gene[i]]
  }
}

# simplify annotations
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
  "Mari_tig00039850-10-191400_plv (rev)",
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
  gene_df_reversed %>% select(molecule, gene, start, end, annotation_short),
  aes(xmin = start, xmax = end, y = molecule, id = annotation_short),
  on = "MCP"
)

ggplot(gene_df_reversed, aes(xmin = start, xmax = end, y = molecule, fill = annotation_short)) +
  geom_gene_arrow(aes(forward = orientation, arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), size = 0.5) +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(
    MCP = "hotpink",            # Stays hotpink - distinct and vibrant
    `DNA Pol` = "#E2D3CEFF",      # Deep, bold red - excellent for highlighted/important
    A32 = "#DFC8CBFF",            # Soft, light blue - subtle, not highlighted
    `GIY-YIG ENase` = "#CD9ABCFF", # Lighter, softer purple
    `His-Me ENase` = "#C28AB1FF",  # Deeper, richer purple - clearly related, but distinct
    `P-loop NTPase` = "#B980A7FF", # Warm, rich burnt orange
    `RNaseH-like sf` = "#AC7299FF", # Deeper, more brownish-orange - clearly related, but distinct
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
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold")
  )

ggsave(plot = last_plot(), file = "final/plv_genome_map.png", height = 4.5, width = 8)
ggsave(plot = last_plot(), file = "final/plv_genome_map.svg", height = 4.5, width = 8)

