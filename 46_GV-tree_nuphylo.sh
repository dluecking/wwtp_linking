#!/bin/bash

# === Step 1: Set up environment ===
echo "🔹 Loading required modules and activating 'nuphylo' Conda environment"
module load Conda
conda activate nuphylo

# Define paths
NCLDV_output_dir="/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking/intermediate/GV_tree/ncldv_output"
NUPHYLO_bin_dir="/lisc/home/user/luecking/bin/NuPhylo"
NUPHYLO_output_dir="/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking/intermediate/GV_tree/nuhpylo_output"

echo "🔹 Creating NuPhylo output directory: $NUPHYLO_output_dir"
mkdir -p "$NUPHYLO_output_dir"

# === Step 2: Link protein directory into NuPhylo directory ===
echo "🔹 Linking NCLDV output directory to NuPhylo bin directory"
ln -sf "$NCLDV_output_dir" "$NUPHYLO_bin_dir/ncldv_output"

# === Step 3: Run NuPhylo ===
echo "🔹 Changing to NuPhylo bin directory: $NUPHYLO_bin_dir"
cd "$NUPHYLO_bin_dir"

# make temporary output directory
mkdir -p tmp_nuphylo_out

# loop over markers, BUT look up the corresponding name
declare -A marker_to_name=(
    [GVOGm0003]="MCP"
    [GVOGm0013]="SF2"
    [GVOGm0022]="RNAPS"
    [GVOGm0023]="RNAPL"
    [GVOGm0054]="PolB"
    [GVOGm0172]="TF2B"
    [GVOGm0461]="Topo2"
    [GVOGm0760]="A32"
    [GVOGm0890]="VLTF3"
)

for marker in GVOGm0003 GVOGm0013 GVOGm0022 GVOGm0023 GVOGm0054 GVOGm0172 GVOGm0461 GVOGm0760 GVOGm0890; do
    echo "🔹 Running NuPhylo for marker: $marker"

    name=${marker_to_name[$marker]}
    python ~/bin/NuPhylo/NuPhylo.py \
        -i ncldv_output/gvs_and_yara_renamed_${marker}.faa \
        -m $name \
        -o tmp_nuphylo_out/${marker}
done

echo "🔹 Moving all NuPhylo output files to final output directory: $NUPHYLO_output_dir"

# cleaning and moving
mv tmp_nuphylo_out/* "$NUPHYLO_output_dir"

unlink ncldv_output
rmdir tmp_nuphylo_out


echo "✅ All NuPhylo calculations completed successfully."
