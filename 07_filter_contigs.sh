#!/bin/bash
module load Conda
conda activate bioinf


# Configuration
INPUT_DIR="intermediate/contigs/lc/"
REMOVE_LIST="helperfiles/contigs_to_remove_since_matching_GVs.txt"

# Ensure seqkit is in PATH (e.g., via active Conda env)
# If not using Conda, ensure seqkit is globally accessible or add its path:
# PATH="/path/to/seqkit/bin:$PATH" # Uncomment and adjust if needed
# conda activate bioinf # Uncomment and adjust if needed

# Validate removal list exists
if [ ! -f "${REMOVE_LIST}" ]; then
    echo "Error: Exclusion list '${REMOVE_LIST}' not found. Please ensure it exists." >&2
    exit 1
fi

# Process each FNA file in the directory
shopt -s nullglob # Handle empty directories gracefully
for FNA_FILE_PATH in "${INPUT_DIR}"*.fna; do
    FILENAME=$(basename "${FNA_FILE_PATH}")
    # Construct output filename with _filtered.fna suffix
    OUTPUT_FILE_PATH="${INPUT_DIR}${FILENAME%.fna}_filtered.fna"

    # Filter out contigs directly using the REMOVE_LIST
    seqkit grep -v -f "${REMOVE_LIST}" "${FNA_FILE_PATH}" > "${OUTPUT_FILE_PATH}" || {
        echo "Error: Filtering failed for ${FNA_FILE_PATH}. Check seqkit output and file permissions." >&2
        exit 1
    }
done
shopt -u nullglob # Turn off nullglob

echo "Filtering complete for contigs in ${INPUT_DIR}."