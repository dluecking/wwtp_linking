# Author: dlu @ veelab
# Version: 2025-05-27

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyr)

# set working directory
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


# summarise ---------------------------------------------------------------

contig_df <- hmm_out %>%
  select(contig, sample, query_type, score) %>%
  # No need for mutate(score = as.numeric(score)) if done above
  group_by(contig, sample, query_type) %>%
  summarise(score = max(score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = query_type, values_from = score, values_fill = NA) %>%
  mutate(
    # Add 'ALL_PLV' to the list of columns for rowSums
    number_of_hits = rowSums(!is.na(across(c(ATPase, MCP, Penton, PRO, ALL_PLV))))
  )

contig_df$contig_ID <- paste0(contig_df$sample, "_", contig_df$contig)


# assign label ------------------------------------------------------------

contig_df <- contig_df %>% 
  filter(number_of_hits >= 3 | !is.na(ALL_PLV))
contig_df$label <- "vph"
contig_df$label[!is.na(contig_df$ALL_PLV)] <- "plv"



# save to helper file -----------------------------------------------------

fwrite(contig_df %>% 
  select(contig, sample, label),
  "helperfiles/vph_and_plvs_to_retrieve.txt")

