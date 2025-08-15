#!/bin/bash
#SBATCH --job-name=prodigal_orf
#SBATCH --output=log/prodigal_output_%A_%a.out
#SBATCH --error=log/prodigal_error_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-62 # 0-62 is the array


module load prodigal


# --- SETUP ---
INPUT_DIR="data/gv_contigs"
FAA_DIR="intermediate/padloc/gv_faa"
GFF_DIR="intermediate/padloc/gv_gff"

# Create an array of the target polished fasta files
FILES=("$INPUT_DIR"/*_polished.fasta)

# Get the specific file for this task using the Slurm array task ID
CURRENT_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

# Define the output file names based on the input file
BASENAME=$(basename "$CURRENT_FILE" .fasta)
OUTPUT_FAA="$FAA_DIR/${BASENAME}.faa"
OUTPUT_GFF="$GFF_DIR/${BASENAME}.gff"

# Ensure output directories2 exist
mkdir -p "$FAA_DIR" "$GFF_DIR"

# --- EXECUTION ---
echo "Running Prodigal on: $CURRENT_FILE"
prodigal -i "$CURRENT_FILE" -p meta -f gff -o "$OUTPUT_GFF" -a "$OUTPUT_FAA"
echo "Done. GFF and FAA files written for $BASENAME"