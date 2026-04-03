# Author: dlu @ veelab
# Version: 2025-06-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)


# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load data and retrieve mapped_to ----------------------------------------

blast_out <- fread("intermediate/integration/blast_results/all_overhangs_blastn.tsv")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                      "qlen", "slen")

blast_out$read <- str_remove(blast_out$qseqid, "\\_to\\_.*$")
blast_out$mapped_to <- str_remove(str_remove(blast_out$qseqid, "^.*\\_to\\_"), "\\_overhang.*$")

# filter out self hits
blast_out <- blast_out %>% 
  filter(mapped_to != sseqid)

blast_out$qcov <- blast_out$length/blast_out$qlen * 100

blast_out8080 <- blast_out %>% filter(pident >= 80, qcov >= 80)


# create two layers for each overhang type --------------------------------

blast_out8080$overhang_type <- if_else(
  stringr::str_ends(blast_out8080$qseqid, "m"), 
  true = "m", 
  false = "b"
)

for(type in unique(blast_out8080$overhang_type)){
  summary_table <- blast_out8080 %>%
    filter(overhang_type == type) %>% 
    group_by(mapped_to, sseqid) %>%
    summarise(number_of_reads = n(), .groups = 'drop') %>%
    rename(from = mapped_to, to = sseqid)
  
  summary_table$log10_reads <- log10(summary_table$number_of_reads + 1)
  names(summary_table) <- c("from", "to", "number_of_reads", paste0("integration_", type))
  # fwrite(summary_table %>% select(from, to, paste0("integration_", type)), paste0("intermediate/network/integration_", type, ".csv"))
}



# can be deleted ----------------------------------------------------------

library(Rsamtools)
library(stringr)
parse_query_coords <- function(cigar, pos, flag) {
  if (is.na(cigar) || is.na(pos) || nchar(cigar) == 0 || cigar == "*") {
    return(list(query_start = NA, query_end = NA, ref_end = NA))
  }
  
  # extract all length-type pairs
  matches <- gregexpr("[0-9]+[A-Z]", cigar)
  tokens  <- regmatches(cigar, matches)[[1]]
  
  if (length(tokens) == 0) return(list(query_start = NA, query_end = NA, ref_end = NA))
  
  lengths <- as.integer(sub("[A-Z]$", "", tokens))
  types   <- sub("^[0-9]+", "", tokens)
  
  query_start <- if (types[1] == "S") lengths[1] + 1 else 1
  query_end   <- sum(lengths[types %in% c("M", "I", "S", "=", "X")])
  ref_end     <- pos + sum(lengths[types %in% c("M", "D", "N", "=", "X")]) - 1
  
  list(query_start = query_start, query_end = query_end, ref_end = ref_end)
}

bam_path <- "intermediate/integration/overhang_mappings/OdNE/overhanging_Aved_tig00303955-10-54120_vph.bam"

param <- ScanBamParam(
  what = c("qname", "flag", "pos", "cigar", "seq"),
  flag = scanBamFlag(isUnmappedQuery = FALSE)
)
# Don't filter df_bam yet — build coords from the FULL bam first
bam <- scanBam(bam_path, param = param)[[1]]

df_bam <- data.frame(
  read    = bam$qname,
  flag    = bam$flag,
  pos     = bam$pos,
  cigar   = bam$cigar,
  seq_len = str_length(bam$seq),
  stringsAsFactors = FALSE
)

df_bam$orientation <- ifelse(bitwAnd(df_bam$flag, 16L) > 0, "reverse", "forward")

# Parse coords on FULL df_bam
coords <- mapply(parse_query_coords,
                 df_bam$cigar, df_bam$pos, df_bam$flag,
                 SIMPLIFY = FALSE)

df_bam$query_start <- sapply(coords, `[[`, "query_start")
df_bam$query_end   <- sapply(coords, `[[`, "query_end")
df_bam$ref_end     <- sapply(coords, `[[`, "ref_end")

# NOW filter
df_bam <- df_bam %>%
  filter(read %in% df$read)

nrow(df_bam)  # should be ~39


# merge both alignments per read into df
a <- df %>%
  select(sseqid, sstart, send, qlen, read_length, read, overhang_type, match_orientation, mapping_orientation) %>%
  left_join(
    df_bam %>% select(read, query_start, query_end, ref_start = pos, ref_end),
    by = "read",
    relationship = "many-to-many"
  ) %>%
  # for each row in df, keep the df_bam alignment whose ref_start is closest to sstart
  mutate(dist = abs(ref_start - sstart)) %>%
  group_by(read, sstart, send) %>%
  slice_min(dist, n = 1) %>%
  ungroup() %>%
  select(-dist)
