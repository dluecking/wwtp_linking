#!/bin/bash
#SBATCH --job-name=cat_and_iqtree
#SBATCH --output=log/REV_cat_and_iqtree_output_%A_%a.out
#SBATCH --error=log/REV_cat_and_iqtree_error_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=36:00:00
#SBATCH --mem=16G


echo "concat"

module load Conda
conda activate catfasta2phyml

NUPHYLO_DIR="intermediate/GV_tree/nuhpylo_output"
OUT_DIR="intermediate/GV_tree/concat_approach"

catfasta2phyml -c -f ${NUPHYLO_DIR}/GVOGm0*/allseqs.trimmed.aln \
  > "${OUT_DIR}/concat.trimmed.fasta" \
  2> "${OUT_DIR}/partitions.txt"

echo "Done with concat"

echo "---"

echo "running iqtree"

module unload Conda
module load IQ-TREE/2.4.0

iqtree2 -T 8 \
  -s "${OUT_DIR}/concat.trimmed.fasta" \
  -B 1000 -alrt 1000 \
  --prefix "${OUT_DIR}/concat_alrt" \
  -redo
  