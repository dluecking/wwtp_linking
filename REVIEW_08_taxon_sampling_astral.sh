#! /bin/bash

module load Conda
conda activate aster-1.23

ASTRAL_dir="intermediate/GV_tree/REVIEW_taxon_sampling/astral"
mkdir -p "$ASTRAL_dir"

# 8 unaffected markers from the REVIEW run (REVIEW_01/02), plus the new
# no-yara GVOGm0760 tree (REVIEW_06/07) in place of the old one
cat intermediate/GV_tree/REVIEW/GVOGm0003.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0013.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0022.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0023.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0054.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0172.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0461.nwk \
    intermediate/GV_tree/REVIEW/GVOGm0890.nwk \
    intermediate/GV_tree/REVIEW_taxon_sampling/iqtree/GVOGm0760.nwk \
    > "$ASTRAL_dir/combined_taxon_sampling.nwk"

astral -i "$ASTRAL_dir/combined_taxon_sampling.nwk" -o "$ASTRAL_dir/combined_taxon_sampling_astral.tre"
