#!/bin/bash
#SBATCH --job-name=run_diamond_tax_array
#SBATCH --output=log/08_diamond_tax_array_%A_%a.out # %A for job ID, %a for array task ID
#SBATCH --error=log/08_diamond_tax_array_%A_%a.err   # Separate logs for each array task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # Match the --threads argument for diamond_assign_classficiation.sh
#SBATCH --time=24:00:00 # Adjust based on expected runtime per file (e.g., 1-12 hours)
#SBATCH --mem=128G # Adjust based on DIAMOND's memory requirements for the NR database
#SBATCH --array=1-1 # up to 23

# load modules
module load DIAMOND/2.1.9-GCC-13.2.0
module load Conda
module load prodigal/2.6.3-GCCcore-13.3.0
conda activate bioinf

# --- Configuration ---
# Path to your diamond_assign_classficiation.sh script
DIAMOND_TAX_SCRIPT="/lisc/home/user/luecking/scripts/diamond_assign_classfication.sh"

# List of input FASTA files (each line should be a full path to a *_filtered.fna file)
INPUT_FILE_LIST="helperfiles/lc_filtered_files.txt"

# Output directory for the taxonomy CSVs
OUTPUT_TAX_DIR="intermediate/lc_tax/"

echo "Starting taxonomic assignment for filtered LC contigs using Slurm Array..."
echo "Input file list: ${INPUT_FILE_LIST}"
echo "Output directory: ${OUTPUT_TAX_DIR}"
echo "Slurm Array Task ID: ${SLURM_ARRAY_TASK_ID}" # Shows which task this specific job is
echo "---"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_TAX_DIR}"

# Validate the diamond_assign_classficiation.sh script exists and is executable
if [ ! -f "${DIAMOND_TAX_SCRIPT}" ]; then
    echo "ERROR: diamond_assign_classficiation.sh script not found at '${DIAMOND_TAX_SCRIPT}'. Exiting." >&2
    exit 1
fi

# Validate the input file list exists
if [ ! -f "${INPUT_FILE_LIST}" ]; then
    echo "ERROR: Input file list '${INPUT_FILE_LIST}' not found. Exiting." >&2
    exit 1
fi

# Determine the input FASTA file for this specific array task
# SLURM_ARRAY_TASK_ID is 0-indexed, so we add 1 for sed's 1-indexing.
INPUT_FASTA_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" "${INPUT_FILE_LIST}")

# Check if INPUT_FASTA_FILE is empty (e.g., if array_task_id exceeds file lines in list)
if [ -z "${INPUT_FASTA_FILE}" ]; then
    echo "WARNING: No input FASTA file found for array task ID ${SLURM_ARRAY_TASK_ID}. This might happen if the array size is larger than the number of files in the list. Exiting this task gracefully." >&2
    exit 0 # Exit successfully for this task, as there's nothing to do
fi

# Extract filename (e.g., "sample_contig_lc_filtered.fna" from "/path/to/sample_contig_lc_filtered.fna")
FILENAME=$(basename "${INPUT_FASTA_FILE}")

# Construct the output CSV filename (e.g., sample_contig_lc_filtered.csv)
OUTPUT_CSV_FILE="${OUTPUT_TAX_DIR}${FILENAME%.fna}.csv"

echo "Processing: ${INPUT_FASTA_FILE}"
echo "Output will be: ${OUTPUT_CSV_FILE}"

# Run the diamond_assign_classficiation.sh script
# It will use the --threads specified by --cpus-per-task in this SBATCH script.
# The assign_taxonomy.sh script handles its own temporary directory creation.
"${DIAMOND_TAX_SCRIPT}" --threads "${SLURM_CPUS_PER_TASK}" \
                        "${INPUT_FASTA_FILE}" "${OUTPUT_CSV_FILE}"

if [ $? -eq 0 ]; then
    echo "Successfully processed ${INPUT_FASTA_FILE} for array task ${SLURM_ARRAY_TASK_ID}."
    rm -rf $TMPDIR
else
    echo "ERROR: Failed to process ${INPUT_FASTA_FILE} for array task ${SLURM_ARRAY_TASK_ID}. Check log for details." >&2
    exit 1 # Exit with error code for this specific array task
fi

echo "Taxonomic assignment for array task ${SLURM_ARRAY_TASK_ID} complete."
