#!/bin/bash
{
  module load Conda
  conda activate bioinf

  for line in $(cat helperfiles/list_of_samples.txt); do
    SHORT_SAMPLE_NAME="${line%%_*}"
    echo "Working on: $SHORT_SAMPLE_NAME"

    python scripts/pysam_filter_and_extract_overhangs.py \
      --input-bam intermediate/integration/combined_mappings/${SHORT_SAMPLE_NAME}_mapped_sorted.bam \
      --fasta-output-dir intermediate/integration/overhang_fasta/${SHORT_SAMPLE_NAME} \
      --bam-output-dir intermediate/integration/overhang_mappings/${SHORT_SAMPLE_NAME}

    echo "Done with: $SHORT_SAMPLE_NAME"
  done
} &> log/25_pysam_filter_extract.out
