#!/bin/bash
#SBATCH --job-name=eggnog
#SBATCH --output=logs/eggnog_output_%A_%a.out
#SBATCH --error=logs/eggnog_error_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=2-63 # up to 63 

# Load necessary modules
module load conda
conda activate eggnog-mapper-2.1.13


# Set variables
FILES_LIST="helperfiles/list_of_gv_protein_files.txt"
INPUT_DIR="intermediate/proteins/gv"
OUTPUT_DIR="intermediate/eggnog"

# Create a temporary directory for each job
mkdir -p $TEMP_DIR

# Get the filename for this specific array task
FILENAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILES_LIST)
BASE_FILENAME=$(basename "$FILENAME" .faa)

echo "Starting emapper.py for file: $FILENAME"
echo "Job ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID"

# Run the emapper.py command
emapper.py -i "$INPUT_DIR/$FILENAME" \
           --output "$OUTPUT_DIR/$BASE_FILENAME" \
           --cpu "$SLURM_CPUS_PER_TASK" \
           --temp_dir=$TMPDIR

# Clean up temporary directory
rm -rf "$TEMP_DIR"

echo "Finished emapper.py for file: $FILENAME"