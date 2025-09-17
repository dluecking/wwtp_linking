#!/bin/bash
#SBATCH --job-name=eggnog
#SBATCH --output=log/eggnog_rerun_output_%A_%a.out
#SBATCH --error=log/eggnog_rerun_error_%A_%a.err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-63 # up to 63

# Load necessary modules
module load conda
conda activate eggnog-mapper-2.1.13


# Set variables
FILES_LIST="helperfiles/list_of_gv_protein_files.txt"
INPUT_DIR="intermediate/proteins/gv"
OUTPUT_DIR="intermediate/eggnog"

# Get the filename for this specific array task
FILENAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILES_LIST)
BASE_FILENAME=$(basename "$FILENAME" .faa)

# Define the output files to check
HITS_FILE="${OUTPUT_DIR}/${BASE_FILENAME}.emapper.hits"
ANNOTATIONS_FILE="${OUTPUT_DIR}/${BASE_FILENAME}.emapper.annotations"

# Check for existing output files and run the mapper if needed
if [[ -f "$HITS_FILE" ]] && [[ -f "$ANNOTATIONS_FILE" ]]; then
    echo "Both files exist for $BASE_FILENAME. Skipping."
elif [[ -f "$HITS_FILE" ]] && [[ ! -f "$ANNOTATIONS_FILE" ]]; then
    echo "Hits file exists but annotations are missing for $BASE_FILENAME. Rerunning."
    emapper.py --override -i "$INPUT_DIR/$FILENAME" --output "$OUTPUT_DIR/$BASE_FILENAME" --cpu "$SLURM_CPUS_PER_TASK" --temp_dir=$TMPDIR 
else
    echo "No files found for $BASE_FILENAME. Running emapper."
    emapper.py -i "$INPUT_DIR/$FILENAME" --output "$OUTPUT_DIR/$BASE_FILENAME" --cpu "$SLURM_CPUS_PER_TASK" --temp_dir=$TMPDIR
fi


echo "Finished emapper.py for file: $FILENAME"