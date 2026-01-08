#!/bin/bash
#SBATCH --job-name=prodigal_predict_4loops
#SBATCH --output=log/prodigal_predict_4loops_%j.out # %j for job ID
#SBATCH --error=log/prodigal_predict_4loops_%j.err   # %j for job ID
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # Prodigal is generally single-threaded per input file
#SBATCH --time=24:00:00 # Adjust based on total processing time for ALL files
#SBATCH --mem=16G # Adjust memory as needed for sequential processing of files

# --- Load Modules ---
module load prodigal/2.6.3-GCCcore-13.3.0

echo "Starting Prodigal protein prediction using four distinct loops..."
echo "---"

# --- Configuration ---
# Base output directory for all predicted proteins
OUTPUT_PROTEIN_BASE_DIR="intermediate/proteins/"

# Ensure the main output directory exists
mkdir -p "${OUTPUT_PROTEIN_BASE_DIR}"

# ==============================================================================
# Loop 1: data/gv_contigs (extension .fasta, output .faa)
# ==============================================================================
CURRENT_INPUT_DIR="data/gv_contigs"
CURRENT_OUTPUT_DIR="${OUTPUT_PROTEIN_BASE_DIR}gv/"

echo "--- Processing 'gv' contigs ---"
echo "Input directory: ${CURRENT_INPUT_DIR}"
echo "Saving predicted proteins to: ${CURRENT_OUTPUT_DIR}"

mkdir -p "${CURRENT_OUTPUT_DIR}"

if [ ! -d "${CURRENT_INPUT_DIR}" ]; then
    echo "ERROR: Input directory '${CURRENT_INPUT_DIR}' not found. Skipping 'gv' processing." >&2
else
    shopt -s nullglob # Enable nullglob to handle cases where no files match
    FASTA_FILES=(${CURRENT_INPUT_DIR}/*.fasta) # Specific to .fasta for gv_contigs
    shopt -u nullglob

    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
        echo "  No .fasta files found in '${CURRENT_INPUT_DIR}'. Skipping 'gv' processing."
    else
        for FASTA_FILE_PATH in "${FASTA_FILES[@]}"; do
            BASENAME=$(basename "$FASTA_FILE_PATH")
            # Remove .fasta extension, then add .faa
            OUT_FILE="${CURRENT_OUTPUT_DIR}${BASENAME%.fasta}.faa"

            echo "  Predicting proteins for: ${FASTA_FILE_PATH} -> ${OUT_FILE}"
            prodigal -i "$FASTA_FILE_PATH" -a "$OUT_FILE" -p meta
            if [ $? -ne 0 ]; then
                echo "  ERROR: Prodigal failed for ${FASTA_FILE_PATH}. Exiting job." >&2
                exit 1
            fi
        done
    fi
fi
echo "Finished 'gv' processing."
echo "---"

# ==============================================================================
# Loop 2: intermediate/contigs/vph (extension .fna, output .faa)
# ==============================================================================
CURRENT_INPUT_DIR="intermediate/contigs/vph"
CURRENT_OUTPUT_DIR="${OUTPUT_PROTEIN_BASE_DIR}vph/"

echo "--- Processing 'vph' contigs ---"
echo "Input directory: ${CURRENT_INPUT_DIR}"
echo "Saving predicted proteins to: ${CURRENT_OUTPUT_DIR}"

mkdir -p "${CURRENT_OUTPUT_DIR}"

if [ ! -d "${CURRENT_INPUT_DIR}" ]; then
    echo "ERROR: Input directory '${CURRENT_INPUT_DIR}' not found. Skipping 'vph' processing." >&2
else
    shopt -s nullglob
    FASTA_FILES=(${CURRENT_INPUT_DIR}/*.fna) # Specific to .fna for vph
    shopt -u nullglob

    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
        echo "  No .fna files found in '${CURRENT_INPUT_DIR}'. Skipping 'vph' processing."
    else
        for FASTA_FILE_PATH in "${FASTA_FILES[@]}"; do
            BASENAME=$(basename "$FASTA_FILE_PATH")
            # Remove .fna extension, then add .faa
            OUT_FILE="${CURRENT_OUTPUT_DIR}${BASENAME%.fna}.faa"

            echo "  Predicting proteins for: ${FASTA_FILE_PATH} -> ${OUT_FILE}"
            prodigal -i "$FASTA_FILE_PATH" -a "$OUT_FILE" -p meta
            if [ $? -ne 0 ]; then
                echo "  ERROR: Prodigal failed for ${FASTA_FILE_PATH}. Exiting job." >&2
                exit 1
            fi
        done
    fi
fi
echo "Finished 'vph' processing."
echo "---"

# ==============================================================================
# Loop 3: intermediate/contigs/lc (extension .fna, output .faa)
# ==============================================================================
CURRENT_INPUT_DIR="intermediate/contigs/lc"
CURRENT_OUTPUT_DIR="${OUTPUT_PROTEIN_BASE_DIR}lc/"

echo "--- Processing 'lc' contigs ---"
echo "Input directory: ${CURRENT_INPUT_DIR}"
echo "Saving predicted proteins to: ${CURRENT_OUTPUT_DIR}"

mkdir -p "${CURRENT_OUTPUT_DIR}"

if [ ! -d "${CURRENT_INPUT_DIR}" ]; then
    echo "ERROR: Input directory '${CURRENT_INPUT_DIR}' not found. Skipping 'lc' processing." >&2
else
    shopt -s nullglob
    FASTA_FILES=(${CURRENT_INPUT_DIR}/*_filtered.fna) # Specific to filtered.fna for lc
    shopt -u nullglob

    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
        echo "  No .fna files found in '${CURRENT_INPUT_DIR}'. Skipping 'lc' processing."
    else
        for FASTA_FILE_PATH in "${FASTA_FILES[@]}"; do
            BASENAME=$(basename "$FASTA_FILE_PATH")
            # Remove .fna extension, then add .faa
            OUT_FILE="${CURRENT_OUTPUT_DIR}${BASENAME%.fna}.faa"

            echo "  Predicting proteins for: ${FASTA_FILE_PATH} -> ${OUT_FILE}"
            prodigal -i "$FASTA_FILE_PATH" -a "$OUT_FILE" -p meta
            if [ $? -ne 0 ]; then
                echo "  ERROR: Prodigal failed for ${FASTA_FILE_PATH}. Exiting job." >&2
                exit 1
            fi
        done
    fi
fi
echo "Finished 'lc' processing."
echo "---"

# ==============================================================================
# Loop 4: intermediate/contigs/plv (extension .fna, output .faa)
# ==============================================================================
CURRENT_INPUT_DIR="intermediate/contigs/plv"
CURRENT_OUTPUT_DIR="${OUTPUT_PROTEIN_BASE_DIR}plv/"

echo "--- Processing 'plv' contigs ---"
echo "Input directory: ${CURRENT_INPUT_DIR}"
echo "Saving predicted proteins to: ${CURRENT_OUTPUT_DIR}"

mkdir -p "${CURRENT_OUTPUT_DIR}"

if [ ! -d "${CURRENT_INPUT_DIR}" ]; then
    echo "ERROR: Input directory '${CURRENT_INPUT_DIR}' not found. Skipping 'plv' processing." >&2
else
    shopt -s nullglob
    FASTA_FILES=(${CURRENT_INPUT_DIR}/*.fna) # Specific to .fna for plv
    shopt -u nullglob

    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
        echo "  No .fna files found in '${CURRENT_INPUT_DIR}'. Skipping 'plv' processing."
    else
        for FASTA_FILE_PATH in "${FASTA_FILES[@]}"; do
            BASENAME=$(basename "$FASTA_FILE_PATH")
            # Remove .fna extension, then add .faa
            OUT_FILE="${CURRENT_OUTPUT_DIR}${BASENAME%.fna}.faa"

            echo "  Predicting proteins for: ${FASTA_FILE_PATH} -> ${OUT_FILE}"
            prodigal -i "$FASTA_FILE_PATH" -a "$OUT_FILE" -p meta
            if [ $? -ne 0 ]; then
                echo "  ERROR: Prodigal failed for ${FASTA_FILE_PATH}. Exiting job." >&2
                exit 1
            fi
        done
    fi
fi
echo "Finished 'plv' processing."
echo "---"

echo "All Prodigal predictions complete."