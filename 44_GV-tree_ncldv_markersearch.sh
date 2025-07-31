#!/bin/bash

# === Step 1: Move GV proteins and Yaravirus proteins to a separate directory ===
echo "🔹 Creating protein input directory: $PROTEIN_DIR"
PROTEIN_DIR="intermediate/GV_tree/ncldv_input_proteins"
mkdir -p $PROTEIN_DIR

echo "🔹 Copying GV and Yaravirus protein files to input directory"
cp intermediate/proteins/gv/*.faa data/MT293574.1.faa $PROTEIN_DIR

# === Step 2: Set up environment ===
echo "🔹 Loading required modules"
module load conda
conda activate ncldv_markersearch
module load hmmer 

# Define paths
NCLDV_bin_dir="/lisc/home/user/luecking/bin/ncldv_markersearch"
NCLDV_output_dir="/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking/intermediate/GV_tree/ncldv_output"

echo "🔹 Creating NCLDV output directory: $NCLDV_output_dir"
mkdir -p $NCLDV_output_dir

# === Step 3: Link protein directory into markersearch directory ===
echo "🔹 Linking protein input directory to ncldv_markersearch bin"
ln -s $PROTEIN_DIR $NCLDV_bin_dir/ncldv_input_proteins

# === Step 4: Run marker search ===
echo "🔹 Changing to ncldv_markersearch bin directory"
cd $NCLDV_bin_dir

for marker in GVOGm0003 GVOGm0013 GVOGm0022 GVOGm0023 GVOGm0054 GVOGm0172 GVOGm0461 GVOGm0760 GVOGm0890; do
    echo "🔹 Running marker search for marker: $marker"

    python ncldv_markersearch.py \
        -i ncldv_input_proteins \
        -n "gvs_and_yara_${marker}" \
        -t 8 \
        -m "$marker" \
        -g muscle

    echo "🔹 Moving output files for $marker to output directory"
    mv gvs_and_yara_${marker}* $NCLDV_output_dir
done

echo "✅ All marker searches completed."