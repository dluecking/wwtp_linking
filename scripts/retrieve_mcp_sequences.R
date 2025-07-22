# Author: dlu @ veelab
# Version: 2025-07-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)


# parsing inputs (gemini) -------------------------------------------------
# Get command-line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)

# Initialize an empty list to store parsed arguments, similar to argparse's 'args' object
args <- list()

i <- 1
while (i <= length(cmd_args)) {
  arg_name <- cmd_args[i]
  if (startsWith(arg_name, "--")) {
    # Remove the "--" prefix to get the clean argument name
    clean_name <- substring(arg_name, 3)
    # Replace hyphens with underscores to match R's typical variable naming
    clean_name <- gsub("-", "_", clean_name)
    
    # Check if there's a corresponding value
    if (i + 1 <= length(cmd_args) && !startsWith(cmd_args[i+1], "--")) {
      args[[clean_name]] <- cmd_args[i+1]
      i <- i + 2 # Move to the next argument pair
    } else {
      # Handle flags that might not have an immediate value (e.g., boolean flags)
      # For this script, all arguments are expected to have a value.
      # If a flag is found without a value, it will be treated as missing.
      warning(paste("Argument", arg_name, "found without a value. Please check your command line."))
      args[[clean_name]] <- NA # Assign NA or NULL if no value is found
      i <- i + 1
    }
  } else {
    # This case should ideally not be hit if arguments are always --flag value
    warning(paste("Unexpected argument format:", arg_name))
    i <- i + 1
  }
}

# Basic validation for required arguments (since commandArgs doesn't do this automatically)
required_args <- c("public_hmmout", "public_proteins", "my_hmmout_dir",
                   "my_vph_protein_dir", "my_plv_protein_dir", "output")

for (req_arg in required_args) {
  if (is.null(args[[req_arg]]) || is.na(args[[req_arg]]) || args[[req_arg]] == "") {
    stop(paste("Error: Required argument --", gsub("_", "-", req_arg), " is missing or empty.", sep=""))
  }
}


# debug mode --------------------------------------------------------------


# Set DEBUG to TRUE to use hardcoded paths for development/testing
DEBUG <- FALSE # Set to TRUE to enable debug mode and overwrite arguments

if (DEBUG) {
  # set working directory
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd("../")
  
  args <- list()
  
  message("DEBUG mode is ON. Overwriting command-line arguments with hardcoded paths.")
  args$public_hmmout       <- "intermediate/hmm_out/all_public_vphs.tbl"
  args$public_proteins     <- "intermediate/proteins/all_public_vphs_proteins.faa"
  args$my_hmmout_dir       <- "intermediate/hmm_out"
  args$my_vph_protein_dir  <- "intermediate/proteins/cleaned/vph"
  args$my_plv_protein_dir  <- "intermediate/proteins/cleaned/plv"
  args$output              <- "intermediate/proteins/all_public_and_my_own_MCP_proteins.faa"
}


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


# load hmm out ------------------------------------------------------------

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
  filter(str_detect(query, "MCP") | query == "ALL_PLV")


# retrieve public MCPs ----------------------------------------------------

all_public_proteins <- read.fasta(args$public_proteins)
all_public_MCPs <- all_public_proteins[unique(hmm_out$target)]

all_public_MCPs_cleaned_1 <- all_public_MCPs[!sapply(all_public_MCPs, is.null)]


# load vph and plv proteins -----------------------------------------------

vph_proteins <- unlist(lapply(list.files(path = args$my_vph_protein_dir, pattern = "\\.faa$", full.names = TRUE), read.fasta), recursive = F)
plv_proteins <- unlist(lapply(list.files(path = args$my_plv_protein_dir, pattern = "\\.faa$", full.names = TRUE), read.fasta), recursive = F)

# find VPH MCPs
vph_MCPs <- list()

for (prot_name in names(vph_proteins)) {
  current_shorted_name <- str_remove(prot_name, "_vph")
  current_shorted_name <- str_remove(current_shorted_name, "^[A-Za-z]*_")
  
  if (current_shorted_name %in% hmm_out$target) {
    vph_MCPs[[prot_name]] <- vph_proteins[[prot_name]]
  }
}

# find PLV MCPs
plv_MCPs <- list()
for (prot_name in names(plv_proteins)) {
  current_shorted_name <- str_remove(prot_name, "_plv")
  current_shorted_name <- str_remove(current_shorted_name, "^[A-Za-z]*_")
  
  if (current_shorted_name %in% hmm_out$target) {
    plv_MCPs[[prot_name]] <- plv_proteins[[prot_name]]
  }
}



# write to file -----------------------------------------------------------

all_MCPs <- c(plv_MCPs, vph_MCPs, all_public_MCPs_cleaned_1)
write.fasta(sequences = all_MCPs,
            names = getName(all_MCPs),
            file.out = args$output)
