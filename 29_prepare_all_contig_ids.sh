#!/bin/bash
# Define input and output paths
INPUT_FASTA="intermediate/contigs/all_contigs.fna"
OUTPUT_DIR="intermediate/network"
OUTPUT_FILE="${OUTPUT_DIR}/all_contig_ids.txt"

mkdir -p "${OUTPUT_DIR}"

# Extract header lines (starting with '>') and remove the leading '>' character
grep '^>' "${INPUT_FASTA}" | sed 's/^>//' > "${OUTPUT_FILE}"

echo "Extracted and cleaned contig IDs from ${INPUT_FASTA} to ${OUTPUT_FILE}"