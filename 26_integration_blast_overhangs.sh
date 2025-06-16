#!/bin/bash
#SBATCH --job-name=integration_blast_overhangs
#SBATCH --output=log/integration_blast_overhangs_%A_%a.out # Standard output log
#SBATCH --error=log/integration_blast_overhangs_%A_%a.err # Standard error log
#SBATCH --ntasks=1 # This script runs sequentially within one task, but can use multiple CPUs
#SBATCH --cpus-per-task=12 # Number of CPUs for BLASTn threads
#SBATCH --time=3-00:00:00  # Max runtime: 3 days
#SBATCH --mem=128G         # Amount of RAM

# --- Load necessary modules ---
module load ncbiblastplus

# --- Define paths ---
DB_DIR="intermediate/contigs"
DB_NAME="all_contigs" # Base name for the BLAST database files
REF_FASTA="${DB_DIR}/${DB_NAME}.fna" # Reference FASTA file for the BLAST database

OVERHANG_FASTA_INPUT_BASE_DIR="intermediate/integration/overhang_fasta" # Base directory for individual overhang FASTA files
ALL_OVERHANGS_FASTA="intermediate/integration/all_overhangs.fasta" # Consolidated FASTA for BLAST query

BLAST_OUTPUT_DIR="intermediate/integration/blast_results" # Directory for BLAST output
BLAST_OUTPUT_FILE="${BLAST_OUTPUT_DIR}/all_overhangs_blastn.tsv" # Final BLAST output file

# Ensure output directories exist
mkdir -p "$DB_DIR"
mkdir -p "$(dirname "$ALL_OVERHANGS_FASTA")"
mkdir -p "$BLAST_OUTPUT_DIR"

# --- Check input reference FASTA existence ---
if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: Reference FASTA file for BLAST database not found: $REF_FASTA" >&2
    exit 1
fi
echo "Reference FASTA file for BLAST database found: $REF_FASTA"

# --- Build BLAST database if it doesn't exist ---
# Check for a common BLAST database file (.nhr) to determine if it's already built
if [ ! -f "${DB_DIR}/${DB_NAME}.nhr" ]; then
    echo "BLAST database for '$DB_NAME' not found. Building..."
    makeblastdb -in "$REF_FASTA" -dbtype nucl -out "${DB_DIR}/${DB_NAME}" -parse_seqids
    echo "BLAST database built successfully."
else
    echo "BLAST database for '$DB_NAME' already exists. Skipping build."
fi

# --- Concatenate all overhang FASTA files into a single query file ---
echo "Concatenating all overhang FASTA files from '$OVERHANG_FASTA_INPUT_BASE_DIR'..."
# Use find to recursively locate all .fasta files and cat them into one
# This assumes files are directly in overhang_fasta or in subdirectories (e.g., sample_ID/file.fasta)
find "$OVERHANG_FASTA_INPUT_BASE_DIR" -type f -name "*.fasta" -print0 | xargs -0 cat > "$ALL_OVERHANGS_FASTA"

# Check if the concatenated FASTA file was created and is not empty
if [ ! -s "$ALL_OVERHANGS_FASTA" ]; then
    echo "ERROR: Concatenated FASTA file '$ALL_OVERHANGS_FASTA' is empty or was not created." >&2
    echo "Please ensure there are .fasta files in '$OVERHANG_FASTA_INPUT_BASE_DIR' or its subdirectories." >&2
    exit 1
fi
echo "All overhang FASTA files concatenated to: $ALL_OVERHANGS_FASTA"

# --- Run BLASTn search ---
echo "Starting BLASTn search..."
blastn -query "$ALL_OVERHANGS_FASTA" \
       -db "${DB_DIR}/${DB_NAME}" \
       -out "$BLAST_OUTPUT_FILE" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
       -num_threads "$SLURM_CPUS_PER_TASK"
echo "BLASTn search completed. Output saved to: $BLAST_OUTPUT_FILE"

echo "Script finished successfully."
