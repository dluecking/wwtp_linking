# this script does only the visualization for 2 selected CRISPR regions and their hits against others
# Author: dlu @ veelab
# Version: 2026-02-27

# Packages
library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)
library(zoo)  # for rollmean

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# functions ---------------------------------------------------------------

#' Smooth coverage with rolling window
#'
#' @param coverage_dt data.table with position and coverage
#' @param window_size Window size for rolling mean
#' @return data.table with smoothed coverage
smooth_coverage <- function(coverage_dt, window_size = 50) {
  coverage_dt[, coverage_smooth := zoo::rollmean(coverage, k = window_size, 
                                                 fill = NA, align = "center")]
  coverage_dt
}


#' Calculate GC content in sliding window
#'
#' @param sequence DNAString or character
#' @param window_size Window size for rolling mean
#' @return data.table with position and GC content
calculate_gc_content <- function(sequence, window_size = 50) {
  seq_char <- as.character(sequence)
  seq_len <- nchar(seq_char)
  
  # Calculate GC at each position (1 = GC, 0 = AT)
  gc_vector <- sapply(1:seq_len, function(i) {
    base <- substr(seq_char, i, i)
    ifelse(base %in% c("G", "C", "g", "c"), 1, 0)
  })
  
  # Rolling mean for smoothing
  gc_rolling <- zoo::rollmean(gc_vector, k = window_size, fill = NA, align = "center")
  
  data.table(
    position = 1:seq_len,
    gc_content = gc_rolling * 100  # Convert to percentage
  )
}

#' Extract coverage from BAM file for a region
#'
#' @param bam_file Path to BAM file
#' @param contig Contig/scaffold name
#' @param start Start position
#' @param end End position
#' @return data.table with position and coverage

extract_coverage <- function(bam_file, contig, start, end) {
  # Create GRanges object for region
  region <- GRanges(seqnames = contig, ranges = IRanges(start = start, end = end))
  
  # Read BAM to get read alignments
  param <- ScanBamParam(which = region, what = c("pos", "qwidth"))
  bam_data <- scanBam(bam_file, param = param)[[1]]
  
  if (length(bam_data$pos) == 0) {
    warning(paste("No reads found in", bam_file, "for", contig, start, "-", end))
    return(data.table(position = start:end, coverage = 0))
  }
  
  # Convert reads to GRanges
  reads_gr <- GRanges(
    seqnames = contig,
    ranges = IRanges(start = bam_data$pos, width = bam_data$qwidth)
  )
  
  # Calculate coverage using GenomicRanges (fast, vectorized)
  cov <- coverage(reads_gr)[[contig]]
  
  # Extract coverage for our region
  if (is.null(cov) || length(cov) < end) {
    warning(paste("Coverage calculation failed for", contig))
    return(data.table(position = start:end, coverage = 0))
  }
  
  coverage_vector <- as.integer(cov[start:end])
  
  data.table(
    position = start:end,
    coverage = coverage_vector
  )
}



# plotting function -------------------------------------------------------

big_plot <- function(fasta_file = fasta_file, ont_bam = ont_bam, ill_bam = ill_bam, contig = contig,
         target_contig = target_contig, start = start, end = end, output_file = output_file,
         gff_file = gff_file, crispr_file = crispr_file, flank_bp = flank_bp, gc_window = gc_window,
         highlight_color = highlight_color, contig_label = contig_label, target_contig_label){
  
  # set color
  HIGHLIGHT_COLOR <- highlight_color
  
  # handle crispr input
  crispr_df <- fread(crispr_file)
  crispr_df$from <- str_remove(crispr_df$`Spacer id`, "\\_CRISPR.*$")
  crispr_df$to <- str_remove(crispr_df$`Target id`, "\\_polypolish$")
  crispr_df <- crispr_df %>% 
    filter(from != to)
  
  # read gff for repeat positions
  gff <- ape::read.gff(gff_file)
  
  # Expand region to include flanks
  region_start <- max(1, start - flank_bp)
  region_end <- end + flank_bp
  
  # Read FASTA sequence
  cat("Reading FASTA...\n")
  fasta <- readDNAStringSet(fasta_file)
  
  # Find matching contig
  contig_idx <- which(names(fasta) == contig)
  if (length(contig_idx) == 0) {
    # Try partial matching
    contig_idx <- grep(contig, names(fasta), fixed = TRUE)
    if (length(contig_idx) == 0) {
      stop(sprintf("Contig '%s' not found in FASTA", contig))
    }
  }
  
  sequence <- fasta[[contig_idx[1]]]
  
  # Check bounds
  if (region_end > length(sequence)) {
    region_end <- length(sequence)
    warning(sprintf("End position adjusted to sequence length: %d", region_end))
  }
  
  # Extract region
  region_seq <- subseq(sequence, region_start, region_end)
  
  # Calculate GC content
  cat("Calculating GC content...\n")
  gc_data <- calculate_gc_content(region_seq, window_size = gc_window)
  gc_data[, position := position + region_start - 1]  # Adjust to genome coordinates
  
  
  # Extract coverage from both BAM files
  cat("Extracting ONT coverage...\n")
  ont_cov <- extract_coverage(ont_bam, contig, region_start, region_end)
  ont_cov <- smooth_coverage(ont_cov, window_size = 50) 
  ont_cov[, read_type := "ONT"]
  
  cat("Extracting Illumina coverage...\n")
  ill_cov <- extract_coverage(ill_bam, contig, region_start, region_end)
  ill_cov <- smooth_coverage(ill_cov, window_size = 50) 
  ill_cov[, read_type := "Illumina"]
  
  # Combine coverage data
  cov_combined <- rbind(ont_cov, ill_cov)
  
  # delete some data
  rm(sequence, ill_cov, ont_cov)
  
  
  # Create data for vertical lines marking ROI boundaries
  roi_boundaries <- data.frame(
    x = c(start, end),
    label = c("ROI Start", "ROI End")
  )
  
  # Highlight region of interest (without flanks) - light grey background
  roi_rect <- data.frame(xmin = start, xmax = end, ymin = -Inf, ymax = Inf)
  
  # filter to only get the correct crispr data
  crispr_df_filtered <- crispr_df %>% filter(to == target_contig & from == contig)
  
  # filter gff
  gff_filtered <- gff %>% filter(seqid == contig & start >= region_start & end <= region_end)
  
  # Keep only repeat_units (drop the repeat_region row)
  repeats <- gff_filtered %>%
    filter(type == "repeat_unit") %>%
    arrange(start) %>%
    mutate(repeat_num = row_number())
  
  # Derive spacers as gaps between consecutive repeats
  spacers <- data.frame(
    start  = repeats$end[-nrow(repeats)] + 1,
    end    = repeats$start[-1] - 1,
    spacer_num = seq_len(nrow(repeats) - 1)
  )
  
  # spacer highlight?
  highlight_spacer <- str_extract(crispr_df_filtered$`Spacer id`, "spacer\\_\\d*$")
  
  # for third plot:
  # Extract spacer number for label
  crispr_df_filtered <- crispr_df_filtered %>%
    mutate(spacer_num_match = as.integer(str_extract(`Spacer id`, "\\d+$")),
           proto_label = paste0("protospacer_", spacer_num_match))
  
  # Rescale protospacer positions to array coordinate space
  array_xmin <- min(repeats$start) - 500
  array_xmax <- max(repeats$end) + 500
  array_range <- array_xmax - array_xmin
  
  proto_scaled <- crispr_df_filtered %>%
    mutate(
      width      = End - Start,
      mid_scaled = scales::rescale((Start + End) / 2,
                                   to = c(array_xmin + array_range * 0.2,
                                          array_xmax - array_range * 0.2)),
      xmin_s = mid_scaled - width / 2,
      xmax_s = mid_scaled + width / 2,
      spacer_num_match = as.integer(str_extract(`Spacer id`, "\\d+$")),
      proto_label = paste0("protospacer_", spacer_num_match)
    )
  
  REPEAT_Y <- 0; SPACER_Y <- 0; BOX_H <- 0.4; PROTO_Y <- -1.2
  
  connector_data <- spacers %>%
    filter(paste0("spacer_", spacer_num) %in% highlight_spacer) %>%
    mutate(x_spacer = (start + end) / 2, y_top = SPACER_Y - BOX_H / 2) %>%
    left_join(proto_scaled %>% select(spacer_num_match, mid_scaled, `N mismatches`),
              by = c("spacer_num" = "spacer_num_match")) %>%
    mutate(x_proto  = mid_scaled,
           y_bottom = PROTO_Y + BOX_H / 2,
           y_mid    = (y_top + y_bottom) / 2,
           mm_label = paste0(`N mismatches`, " MM"))
  
  sz   <- 10 / 2.845  # ~ 3.52 — all non-title text
  sz_t <- 12 / 2.845  # ~ 4.22 — unused here but consistent
  
  # ── p1: GC content ───────────────────────────────────────────────────────────
  p1 <- ggplot(gc_data, aes(x = position, y = gc_content)) +
    geom_rect(data = roi_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey90", alpha = 0.8, inherit.aes = FALSE) +
    geom_vline(data = roi_boundaries, aes(xintercept = x),
               linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_line(color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 50, linetype = "dotted", color = "gray50", alpha = 0.5) +
    # annotate("text", x = start, y = max(gc_data$gc_content, na.rm = TRUE) * 0.95,
    #          label = "Start", hjust = 1.1, vjust = 1, size = sz, color = "red") +
    # annotate("text", x = end,   y = max(gc_data$gc_content, na.rm = TRUE) * 0.95,
    #          label = "End",   hjust = -0.1, vjust = 1, size = sz, color = "red") +
    scale_x_continuous(limits = c(region_start, region_end)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = paste0(contig_label, " targets ", target_contig_label), x = NULL, y = "GC (%)") +
    theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold"))
  
  # ── p2: Coverage ─────────────────────────────────────────────────────────────
  p2 <- ggplot(cov_combined, aes(x = position, y = coverage_smooth, linetype = read_type)) +
    geom_rect(data = roi_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey90", alpha = 0.8, inherit.aes = FALSE) +
    geom_vline(data = roi_boundaries, aes(xintercept = x),
               linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_line(color = "black", linewidth = 0.8) +
    scale_linetype_manual(values = c("ONT" = "dashed", "Illumina" = "dotted")) +
    scale_x_continuous(limits = c(region_start, region_end),
                       labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Genomic Position (bp)", y = "Coverage") +
    theme_cowplot(font_size = 10) +
    theme(legend.position = "none")
  
  # ── p3: CRISPR array + protospacer track ─────────────────────────────────────
  p3 <- ggplot() +
    geom_segment(aes(x = array_xmin, xend = array_xmax, y = 0,       yend = 0),       color = "black", linewidth = 1.5) +
    geom_segment(aes(x = array_xmin, xend = array_xmax, y = PROTO_Y, yend = PROTO_Y), color = "black", linewidth = 1.5) +
    geom_rect(data = repeats,
              aes(xmin = start, xmax = end, ymin = REPEAT_Y - BOX_H/2, ymax = REPEAT_Y + BOX_H/2),
              fill = "grey40", color = "black", linewidth = 0.4) +
    geom_rect(data = spacers,
              aes(xmin = start, xmax = end, ymin = SPACER_Y - BOX_H/2, ymax = SPACER_Y + BOX_H/2,
                  fill = paste0("spacer_", spacer_num) %in% highlight_spacer),
              color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = HIGHLIGHT_COLOR), guide = "none") +
    geom_text(data = spacers %>% filter(paste0("spacer_", spacer_num) %in% highlight_spacer),
              aes(x = (start + end) / 2, y = SPACER_Y + BOX_H/2 + 0.1,
                  label = paste0("spacer_", spacer_num)),
              size = sz, vjust = 0, color = HIGHLIGHT_COLOR, fontface = "bold") +
    # annotate("text", x = min(repeats$start), y = SPACER_Y + BOX_H/2 + 0.15,
    #          label = scales::comma(start), size = sz, hjust = 1, color = "black") +
    # annotate("text", x = max(repeats$end), y = SPACER_Y + BOX_H/2 + 0.15,
    #          label = scales::comma(end),   size = sz, hjust = 0, color = "black") +
    geom_rect(data = proto_scaled,
              aes(xmin = xmin_s, xmax = xmax_s, ymin = PROTO_Y - BOX_H/2, ymax = PROTO_Y + BOX_H/2),
              fill = HIGHLIGHT_COLOR, color = "black", linewidth = 0.4) +
    geom_text(data = proto_scaled,
              aes(x = mid_scaled, y = PROTO_Y - BOX_H/2 - 0.08, label = proto_label),
              size = sz, vjust = 1, color = HIGHLIGHT_COLOR, fontface = "bold") +
    geom_segment(data = connector_data,
                 aes(x = x_spacer, xend = x_proto, y = y_top, yend = y_bottom),
                 color = HIGHLIGHT_COLOR, linewidth = 0.6, linetype = "dashed") +
    geom_label(data = connector_data,
               aes(x = (x_spacer + x_proto) / 2, y = y_mid, label = mm_label),
               size = sz, color = HIGHLIGHT_COLOR,
               fill = "white", label.size = 0.25, label.padding = unit(0.15, "lines")) +
    annotate("text", x = array_xmin - 100, y = BOX_H/2 + 0.3,            label = contig_label,
             size = sz, hjust = 0, color = "black", fontface = "bold") +
    annotate("text", x = array_xmin - 100, y = PROTO_Y + BOX_H/2 + 0.3,  label = target_contig_label,
             size = sz, hjust = 0, color = "black", fontface = "bold") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    coord_cartesian(xlim = c(array_xmin, array_xmax), ylim = c(-2.0, 0.8))
  
  plot_final <- (p1 / p2 / p3) + 
    plot_layout(heights = c(1, 1, 2))
  
  
  ggsave(plot = plot_final, file = output_file, height = 5, width = 7)
  ggsave(plot = plot_final, file = str_replace(output_file, "\\.png", "\\.svg"), height = 5, width = 7)
  ggsave(plot = plot_final, file = str_replace(output_file, "\\.png", "\\.pdf"), height = 5, width = 7)
}



# plot --------------------------------------------------------------------

# arguments, but manually
fasta_file <- "data/gv_contigs/Lyne_1_polished.fasta"
ont_bam <- "intermediate/occurance/mapping/ont/SRR11673972_mapped_sorted.bam"
ill_bam <- "intermediate/occurance/mapping/illumina/SRR11674015_mapped_sorted.bam"
contig <- "Lyne_1"
target_contig <- "Rand_1_2"
start <- 422822
end <- 424179
output_file <- "final/repeat_region_plots/Lyne_1_to_Rand_1_2_gc_cov_array.png"
gff_file <- "intermediate/minced/minced_results.gff"
crispr_file <- "intermediate/CRISPR/map_spacers_to_targets/minced_results_spacers_filtered_lc_vs_target_db_all_hits.tsv"
flank_bp <- 500
gc_window <- 50
highlight_color <- "steelblue"
contig_label <- "Lyne_NCV-Im_1"
target_contig_label <- "Rand_NCV-Im_1"

big_plot(fasta_file <- fasta_file, ont_bam = ont_bam, ill_bam = ill_bam, contig = contig,
         target_contig = target_contig, start = start, end = end, output_file = output_file,
         gff_file = gff_file, crispr_file = crispr_file, flank_bp = flank_bp, gc_window = gc_window, 
         highlight_color = highlight_color, contig_label = contig_label, target_contig_label = target_contig_label)

# arguments, but manually
fasta_file <- "intermediate/contigs/lc/Hade_lc_filtered.fna"
ont_bam <- "intermediate/occurance/mapping/ont/SRR11673977_mapped_sorted.bam"
ill_bam <- "intermediate/occurance/mapping/illumina/SRR11674031_mapped_sorted.bam"
contig <- "Hade_tig01683932-10-1993460_lc"
target_contig <- "Hjor_1"
start <- 132157
end <- 134545
output_file <- "final/repeat_region_plots/Hade_tig01683932-10-1993460_lc_to_Hjor_1_gc_cov_array.png"
gff_file <- "intermediate/minced/minced_results.gff"
crispr_file <- "intermediate/CRISPR/map_spacers_to_targets/minced_results_spacers_filtered_lc_vs_target_db_all_hits.tsv"
flank_bp <- 500
gc_window <- 50
highlight_color <- "seagreen"
contig_label <- "Hade_tig01683932"
target_contig_label <- "Hjor_NCV-Pi_1"

big_plot(fasta_file <- fasta_file, ont_bam = ont_bam, ill_bam = ill_bam, contig = contig,
         target_contig = target_contig, start = start, end = end, output_file = output_file,
         gff_file = gff_file, crispr_file = crispr_file, flank_bp = flank_bp, gc_window = gc_window, 
         highlight_color = highlight_color, contig_label = contig_label, target_contig_label = target_contig_label)


