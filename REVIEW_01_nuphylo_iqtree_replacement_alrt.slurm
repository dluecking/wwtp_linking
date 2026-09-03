#!/bin/bash
#SBATCH --job-name=iqtree2
#SBATCH --output=log/REV_alrt_iqtree2_output_%A_%a.out
#SBATCH --error=log/REV_alrt_iqtree2_error_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --array=0-8

module load IQ-TREE/2.4.0

# same 9 GVOG markers used by 46_GV-tree_nuphylo.sh
markers=(GVOGm0003 GVOGm0013 GVOGm0022 GVOGm0023 GVOGm0054 GVOGm0172 GVOGm0461 GVOGm0760 GVOGm0890)
marker=${markers[$SLURM_ARRAY_TASK_ID]}

NUPHYLO_output_dir="intermediate/GV_tree/nuhpylo_output"
REVIEW_output_dir="intermediate/GV_tree/REVIEW_alrt"

mkdir -p "$REVIEW_output_dir"

aln="${NUPHYLO_output_dir}/${marker}/allseqs.trimmed.aln"

echo "Running IQ-TREE2 for marker: $marker"
iqtree2 -T 8 -s "$aln" -B 1000 -alrt 1000 --prefix "${REVIEW_output_dir}/${marker}_alrt" -redo

# keep a plain .nwk copy alongside the .treefile so the astral step
# (REVIEW_02_astral_on_new_trees.sh) can just glob *.nwk
cp "${REVIEW_output_dir}/${marker}_alrt.treefile" "${REVIEW_output_dir}/${marker}_alrt.nwk"

echo "IQ-TREE2 finished for marker: $marker"

