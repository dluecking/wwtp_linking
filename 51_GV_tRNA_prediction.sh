#!/bin/bash
#SBATCH --job-name=aragorn_tRNA
#SBATCH --output=log/aragorn_output_%A_%a.out
#SBATCH --error=log/aragorn_error_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --array=1-62 # 0 to 62 is the array

module load aragorn-1.2.41

# --- SETUP ---
INPUT_DIR="data/gv_contigs"
OUTPUT_DIR="intermediate/tRNA_prediction"

# Create an array of the target polished fasta files
FILES=("$INPUT_DIR"/*_polished.fasta)

# Get the specific file for this task using the Slurm array task ID
CURRENT_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

# Define the output file name based on the input file
BASENAME=$(basename "$CURRENT_FILE" _polished.fasta)
OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}.trnas.txt"


# --- EXECUTION ---
echo "Running ARAGORN on: $CURRENT_FILE"
aragorn -t -i3000 -w -o "$OUTPUT_FILE" "$CURRENT_FILE"
echo "Done. Output written to: $OUTPUT_FILE"