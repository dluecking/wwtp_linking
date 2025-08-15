#!/bin/bash
#SBATCH --job-name=padloc
#SBATCH --output=log/padloc_output_%A_%a.out
#SBATCH --error=log/padloc_error_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-62%15 # To get this number, run: ls -1 data/gv_contigs/*_polished.fasta | wc -l

module load conda
conda activate padloc


# --- SETUP ---
# Define directories
CONTIG_DIR="data/gv_contigs"
FAA_DIR="intermediate/padloc/gv_faa"
GFF_DIR="intermediate/padloc/gv_gff"
OUT_DIR="intermediate/padloc/padloc_out"

# Create an array of the original polished fasta files to easily get basenames
FILES=("$CONTIG_DIR"/*_polished.fasta)

# Get the basename for this task
CURRENT_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$CURRENT_FILE" .fasta)

# Define the specific input files and the output directory for this task
CURRENT_FAA="$FAA_DIR/${BASENAME}.faa"
CURRENT_GFF="$GFF_DIR/${BASENAME}.gff"
CURRENT_OUT_DIR="$OUT_DIR/$BASENAME"

# Ensure output directories exist
mkdir -p "$CURRENT_OUT_DIR"

# --- EXECUTION ---
echo "Running PADLOC on: $BASENAME"
padloc --faa "$CURRENT_FAA" --gff "$CURRENT_GFF" --outdir "$CURRENT_OUT_DIR" --fix-prodigal
echo "Done. PADLOC results written to: $CURRENT_OUT_DIR"