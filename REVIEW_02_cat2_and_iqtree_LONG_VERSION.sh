#!/bin/bash
#SBATCH --job-name=LV_cat_and_iqtree
#SBATCH --output=log/REV_LV_cat_and_iqtree_output_%A_%a.out
#SBATCH --error=log/REV_LV_cat_and_iqtree_error_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=32G


echo "concat"

module load Conda
conda activate catfasta2phyml

NUPHYLO_DIR="intermediate/GV_tree/nuhpylo_output"
OUT_DIR="intermediate/GV_tree/concat_approach"

catfasta2phyml -c -f ${NUPHYLO_DIR}/GVOGm0*/allseqs.trimmed.aln \
  > "${OUT_DIR}/LONG_VERSION_concat.trimmed.fasta" \
  2> "${OUT_DIR}/LONG_VERSION_partitions.txt"

echo "Done with concat"

echo "---"

echo "running iqtree"

module unload Conda
module load IQ-TREE/2.4.0

iqtree2 -T 8 \
  -s "${OUT_DIR}/LONG_VERSION_concat.trimmed.fasta" \
  -b 1000 \
  --prefix "${OUT_DIR}/LONG_VERSION_concat_alrt" \
  -redo
  