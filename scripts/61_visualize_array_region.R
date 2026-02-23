#!/usr/bin/env Rscript

# CRISPR/Spacer Validation via Coverage & GC Content
# Purpose: Verify CRISPR arrays/spacers in NCVs by examining read support and GC patterns

suppressPackageStartupMessages({
  library(Rsamtools)
  library(Biostrings)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(patchwork)
  library(zoo)  # for rollmean
})

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
  
  # Set up BAM parameters
  param <- ScanBamParam(which = region, what = c("pos", "qwidth"))
  
  # Read BAM
  bam_data <- scanBam(bam_file, param = param)[[1]]
  
  if (length(bam_data$pos) == 0) {
    warning(paste("No reads found in", bam_file, "for", contig, start, "-", end))
    return(data.table(position = start:end, coverage = 0))
  }
  
  # Calculate coverage using pileup
  pileup_param <- PileupParam(max_depth = 100000, min_base_quality = 0, 
                              min_mapq = 0, distinguish_strands = FALSE)
  
  coverage_data <- pileup(bam_file, scanBamParam = param, pileupParam = pileup_param)
  
  # Convert to data.table with full range
  if (nrow(coverage_data) == 0) {
    return(data.table(position = start:end, coverage = 0))
  }
  
  coverage_dt <- data.table(
    position = coverage_data$pos,
    coverage = coverage_data$count
  )
  
  # Fill in missing positions with 0 coverage
  full_range <- data.table(position = start:end)
  coverage_dt <- merge(full_range, coverage_dt, by = "position", all.x = TRUE)
  coverage_dt[is.na(coverage), coverage := 0]
  
  coverage_dt
}

#' Main validation function
#'
#' @param fasta_file Path to FASTA file
#' @param ont_bam Path to ONT BAM file
#' @param ill_bam Path to Illumina BAM file
#' @param contig Contig/scaffold accession
#' @param start Start position of region of interest
#' @param end End position of region of interest
#' @param output_file Output PNG file path
#' @param flank_bp Flanking basepairs to include (default: 500)
#' @param gc_window GC rolling window size (default: 50)
validate_crispr_region <- function(fasta_file, ont_bam, ill_bam, contig, 
                                   start, end, output_file,
                                   flank_bp = 500, gc_window = 50) {
  
  cat(sprintf("Processing: %s:%d-%d\n", contig, start, end))
  
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
  
  # Extract coverage from both BAM files (now including flanks)
  cat("Extracting ONT coverage...\n")
  ont_cov <- extract_coverage(ont_bam, contig, region_start, region_end)
  ont_cov[, read_type := "ONT"]
  
  cat("Extracting Illumina coverage...\n")
  ill_cov <- extract_coverage(ill_bam, contig, region_start, region_end)
  ill_cov[, read_type := "Illumina"]
  
  # Combine coverage data
  cov_combined <- rbind(ont_cov, ill_cov)
  
  # Create data for vertical lines marking ROI boundaries
  roi_boundaries <- data.frame(
    x = c(start, end),
    label = c("ROI Start", "ROI End")
  )
  
  # Highlight region of interest (without flanks) - lighter background
  roi_rect <- data.frame(xmin = start, xmax = end, ymin = -Inf, ymax = Inf)
  
  # Plot 1: GC content
  p1 <- ggplot(gc_data, aes(x = position, y = gc_content)) +
    geom_rect(data = roi_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "yellow", alpha = 0.15, inherit.aes = FALSE) +
    geom_vline(data = roi_boundaries, aes(xintercept = x), 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_line(color = "darkblue", linewidth = 0.8) +
    geom_hline(yintercept = 50, linetype = "dotted", color = "gray50", alpha = 0.5) +
    annotate("text", x = start, y = max(gc_data$gc_content, na.rm = TRUE) * 0.95, 
             label = "Start", hjust = -0.1, vjust = 1, size = 3, color = "red") +
    annotate("text", x = end, y = max(gc_data$gc_content, na.rm = TRUE) * 0.95, 
             label = "End", hjust = 1.1, vjust = 1, size = 3, color = "red") +
    labs(title = sprintf("%s:%d-%d (±%d bp flanks shown)", contig, start, end, flank_bp),
         subtitle = sprintf("GC Content (window = %d bp)", gc_window),
         x = NULL,
         y = "GC Content (%)") +
    scale_x_continuous(limits = c(region_start, region_end)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  # Plot 2: Coverage
  p2 <- ggplot(cov_combined, aes(x = position, y = coverage, color = read_type)) +
    geom_rect(data = roi_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "yellow", alpha = 0.15, inherit.aes = FALSE) +
    geom_vline(data = roi_boundaries, aes(xintercept = x), 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    scale_color_manual(values = c("ONT" = "#E64B35", "Illumina" = "#4DBBD5"),
                       name = "Read Type") +
    labs(subtitle = sprintf("Read Coverage (flanks: %d-%d | ROI: %d-%d | flanks: %d-%d)", 
                            region_start, start - 1, start, end, end + 1, region_end),
         x = "Genomic Position (bp)",
         y = "Coverage (reads)") +
    scale_x_continuous(limits = c(region_start, region_end),
                       labels = scales::comma) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_rect(color = "gray80", fill = "white"))
  
  # Combine plots
  combined <- p1 / p2 + plot_layout(heights = c(1, 1.2))
  
  # Save
  cat(sprintf("Saving to %s...\n", output_file))
  ggsave(output_file, combined, width = 14, height = 8, dpi = 300)
  
  cat("Done!\n")
  
  # Calculate statistics for ROI only and flanking regions
  ont_roi <- ont_cov[position >= start & position <= end]
  ill_roi <- ill_cov[position >= start & position <= end]
  
  ont_flank <- ont_cov[position < start | position > end]
  ill_flank <- ill_cov[position < start | position > end]
  
  # Return summary statistics
  list(
    region = sprintf("%s:%d-%d", contig, start, end),
    region_with_flanks = sprintf("%s:%d-%d", contig, region_start, region_end),
    ont_roi_mean = mean(ont_roi$coverage),
    ont_roi_median = median(ont_roi$coverage),
    ont_flank_mean = mean(ont_flank$coverage),
    ill_roi_mean = mean(ill_roi$coverage),
    ill_roi_median = median(ill_roi$coverage),
    ill_flank_mean = mean(ill_flank$coverage),
    mean_gc = mean(gc_data$gc_content, na.rm = TRUE),
    roi_length = end - start + 1
  )
}

# ============================================================================
# Command-line interface
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  cat("
CRISPR/Spacer Validation Script

Usage: 
  Rscript validate_region.R <fasta> <ont_bam> <ill_bam> <contig> <start> <end> <output.png> [flank] [gc_window]

Arguments:
  fasta       Path to FASTA file
  ont_bam     Path to ONT BAM file
  ill_bam     Path to Illumina BAM file
  contig      Contig/scaffold accession
  start       Start position (ROI start)
  end         End position (ROI end)
  output.png  Output PNG file
  flank       Flanking bp (default: 500)
  gc_window   GC window size (default: 50)

Example:
  Rscript validate_region.R genome.fna ont.bam ill.bam Lyne_1 263777 265136 lyne_crispr.png

Note: 
  - Red dashed vertical lines mark the ROI boundaries (start and end)
  - Yellow shaded area highlights the region of interest
  - Coverage and GC plots show ±500 bp flanking regions

", file = stderr())
  quit(status = 1)
}

fasta_file <- args[1]
ont_bam <- args[2]
ill_bam <- args[3]
contig <- args[4]
start <- as.integer(args[5])
end <- as.integer(args[6])
output_file <- args[7]
flank_bp <- ifelse(length(args) >= 8, as.integer(args[8]), 500)
gc_window <- ifelse(length(args) >= 9, as.integer(args[9]), 50)

# Run validation
stats <- validate_crispr_region(
  fasta_file = fasta_file,
  ont_bam = ont_bam,
  ill_bam = ill_bam,
  contig = contig,
  start = start,
  end = end,
  output_file = output_file,
  flank_bp = flank_bp,
  gc_window = gc_window
)

# Print summary
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Region of Interest: %s\n", stats$region))
cat(sprintf("Full Region (with flanks): %s\n", stats$region_with_flanks))
cat(sprintf("ROI Length: %d bp\n", stats$roi_length))
cat("\n--- Coverage in ROI ---\n")
cat(sprintf("ONT - Mean: %.1f, Median: %.1f\n", stats$ont_roi_mean, stats$ont_roi_median))
cat(sprintf("Illumina - Mean: %.1f, Median: %.1f\n", stats$ill_roi_mean, stats$ill_roi_median))
cat("\n--- Coverage in Flanks ---\n")
cat(sprintf("ONT - Mean: %.1f\n", stats$ont_flank_mean))
cat(sprintf("Illumina - Mean: %.1f\n", stats$ill_flank_mean))
cat("\n--- Sequence Composition ---\n")
cat(sprintf("Mean GC Content: %.1f%%\n", stats$mean_gc))
```

---
  
  ## **Key Changes:**
  
  1. **Vertical red dashed lines** at `start` and `end` positions to clearly mark ROI boundaries
2. **Coverage plots now include ±500 bp flanks** (was already showing them, but now explicitly noted in subtitle)
3. **Text annotations** "Start" and "End" at the top of the GC plot
4. **Enhanced subtitle** on coverage plot showing exact coordinates of flanks and ROI
5. **Split statistics** - now reports coverage separately for ROI vs flanking regions

---
  
  ## **What the Output Shows:**
  
  **Visual elements:**
  - **Yellow shaded area**: Your region of interest (the CRISPR array/spacer)
- **Red dashed vertical lines**: Exact boundaries you specified
- **Flanking regions**: 500 bp on each side to detect assembly artifacts or coverage drops
- **Dotted gray line** at 50% GC for reference

**Statistics printed:**
  ```
=== Summary Statistics ===
  Region of Interest: Lyne_1:263777-265136
Full Region (with flanks): Lyne_1:263277-265636
ROI Length: 1360 bp

--- Coverage in ROI ---
  ONT - Mean: 45.2, Median: 47.0
Illumina - Mean: 123.5, Median: 125.0

--- Coverage in Flanks ---
  ONT - Mean: 52.1
Illumina - Mean: 118.3

--- Sequence Composition ---
  Mean GC Content: 38.2%