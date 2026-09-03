#!/bin/bash

# === Set up environment ===
module load Conda
conda activate nuphylo

NUPHYLO_bin_dir="/lisc/home/user/luecking/bin/NuPhylo"
REVIEW_dir="/lisc/home/user/luecking/luecking_scratch/projects/wwtp_linking/intermediate/GV_tree/REVIEW_taxon_sampling/nuphylo"

# only GVOGm0760 (A32) is affected by yara removal
marker="GVOGm0760"

ln -sf "$REVIEW_dir" "$NUPHYLO_bin_dir/review_taxon_sampling"

cd "$NUPHYLO_bin_dir"
mkdir -p tmp_nuphylo_out

python ~/bin/NuPhylo/NuPhylo.py \
    -i review_taxon_sampling/gvs_and_yara_removed_${marker}.faa \
    -m A32 \
    -o tmp_nuphylo_out/${marker}

mv tmp_nuphylo_out/${marker} "$REVIEW_dir/"

unlink review_taxon_sampling
rmdir tmp_nuphylo_out
