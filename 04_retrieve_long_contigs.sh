#!/bin/bash
#SBATCH --job-name=filter_contigs_by_len
#SBATCH --output=log/04_filter_contigs_by_len_%A.out
#SBATCH --error=log/04_filter_contigs_by_len_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=8G

# --- Configuration ---
ASSEMBLY_DIR="data/assemblies/"
OUTPUT_DIR="intermediate/contigs/lc/"
MIN_LENGTH=100000

# --- Load tools ---
module load Conda
conda activate bioinf # Ensure seqkit is in this environment

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Starting contig filtering and renaming..."

# Loop through each assembly file
for ASSEMBLY_FILE_PATH in "${ASSEMBLY_DIR}"/*.fa; do
    SAMPLE_NAME=$(basename "${ASSEMBLY_FILE_PATH}" | cut -d'_' -f1)
    OUTPUT_FILE_PATH="${OUTPUT_DIR}${SAMPLE_NAME}_lc.fna" # Changed to single output file per sample

    echo "Processing ${ASSEMBLY_FILE_PATH}"

    # Use awk to filter by length and rename headers, then pipe to seqkit for wrapping
    awk -v min_len="${MIN_LENGTH}" \
        -v sample_name="${SAMPLE_NAME}" \
        '
        BEGIN {
            RS=">"; FS="\n"
            OFS="\n"
        }
        NR > 1 {
            header_line = $1
            match(header_line, /^([^ ]+)/, id_match)
            contig_id = id_match[1]

            sequence = ""
            for (i = 2; i <= NF; i++) {
                sequence = sequence $i
            }
            gsub(/[^A-Za-z]/, "", sequence) # Remove non-alphabet characters (including newlines)

            if (length(sequence) >= min_len) {
                new_accession = sample_name "_" contig_id "_lc"
                print ">" new_accession
                print sequence
            }
        }
    ' "${ASSEMBLY_FILE_PATH}" | seqkit seq -w 60 > "${OUTPUT_FILE_PATH}"

    if [ $? -eq 0 ]; then
        echo "  Saved filtered and renamed contigs to: ${OUTPUT_FILE_PATH}"
    else
        echo "  ERROR processing ${ASSEMBLY_FILE_PATH}"
    fi
done

echo "Contig filtering and renaming complete."
