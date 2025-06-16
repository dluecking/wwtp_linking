#!/bin/bash
#SBATCH --job-name=retrieve_contigs
#SBATCH --output=log/retrieve_contigs_%A.out
#SBATCH --error=log/retrieve_contigs_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # seqkit grep is usually fast, 1 CPU is often enough. Adjust if needed.
#SBATCH --time=1:00:00 # Adjust time based on the number/size of contigs
#SBATCH --mem=4G # Adjust memory based on assembly file sizes

# --- Configuration ---
HELPER_FILE="helperfiles/vph_and_plvs_to_retrieve.txt"
ASSEMBLY_DIR="data/assemblies/"
BASE_OUTPUT_DIR="intermediate/contigs/" # Base directory for all retrieved contigs

# --- Load necessary tools ---
module load conda        # Load the conda module if your system uses it
conda activate bioinf    # Activate your specific bioinformatics Conda environment (e.g., where seqkit is installed)

echo "Starting contig retrieval process..."
echo "Helper file: ${HELPER_FILE}"
echo "Assembly directory: ${ASSEMBLY_DIR}"
echo "Base output directory: ${BASE_OUTPUT_DIR}"
echo "---"

# Create the base output directory if it doesn't exist.
# Specific 'vph' and 'plv' subdirectories will be created on the fly.
mkdir -p "${BASE_OUTPUT_DIR}"

# Skip the header line in the helper file and process each entry
tail -n +2 "${HELPER_FILE}" | while IFS=',' read -r contig_name sample_name label; do
    echo "Processing: Contig=${contig_name}, Sample=${sample_name}, Label=${label}"

    # Determine the specific output subdirectory based on the label (vph or plv)
    # This also handles cases where label might not be exactly 'vph' or 'plv' (though your sample data suggests it will be)
    CURRENT_OUTPUT_DIR="${BASE_OUTPUT_DIR}${label}/"

    # Create the specific output subdirectory if it doesn't exist
    mkdir -p "${CURRENT_OUTPUT_DIR}"

    # Construct the expected assembly file name using 'find'
    ASSEMBLY_FILE=$(find "${ASSEMBLY_DIR}" -maxdepth 1 -name "${sample_name}_*assembly_noNs_renamed.fa" | head -n 1)

    if [ -z "${ASSEMBLY_FILE}" ]; then
        echo "  ERROR: Assembly file for sample '${sample_name}' not found in ${ASSEMBLY_DIR}. Skipping."
        continue
    fi

    echo "  Found assembly: ${ASSEMBLY_FILE}"

    # Define the output file name for the extracted contig
    OUTPUT_CONTIG_PATH="${CURRENT_OUTPUT_DIR}${sample_name}_${contig_name}_${label}.fna"
    
    # Define the new accession for the contig within the output file
    NEW_ACCESSION="${sample_name}_${contig_name}_${label}"

    # Use seqkit grep to extract the contig and then sed to rename its header
    seqkit grep -p "${contig_name}" "${ASSEMBLY_FILE}" | \
    sed "s/^>.*/>${NEW_ACCESSION}/" > "${OUTPUT_CONTIG_PATH}"

    if [ $? -eq 0 ]; then
        echo "  Successfully extracted and saved to: ${OUTPUT_CONTIG_PATH}"
    else
        echo "  ERROR: Failed to extract contig ${contig_name} from ${ASSEMBLY_FILE}. Check seqkit and sed output."
    fi
    echo "---"
done

echo "Contig retrieval complete."