#! /bin/bash

module load Conda
conda activate aster-1.23

# concat the REVIEW iqtree gene trees (one .nwk per marker, from REVIEW_01_nuphylo_iqtree_replacement.slurm)
cat intermediate/GV_tree/REVIEW_alrt/*alrt.nwk > intermediate/GV_tree/REVIEW_alrt/combined_review_alrt.nwk

# run the astral command
astral -i intermediate/GV_tree/REVIEW_alrt/combined_review_alrt.nwk -o intermediate/GV_tree/REVIEW_alrt/combined_review_alrt_astral.nwk
