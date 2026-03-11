#!/bin/bash
module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R

# 61_viz_repeat_regions.sh
# Visualize CRISPR repeat regions with coverage and GC content

set -e  # Exit on error

# Define base directories
FASTA_DIR="data/gv_contigs"
ONT_DIR="intermediate/occurance/mapping/ont"
ILL_DIR="intermediate/occurance/mapping/illumina"
OUT_DIR="final/repeat_region_plots"

# Create output directory
mkdir -p ${OUT_DIR}

# SRR mapping: Station -> Illumina,ONT
declare -A SRR_MAP=(
    ["EsbE"]="SRR11674041,SRR11673980"
    ["Hjor"]="SRR11674022,SRR11673974"
    ["Lyne"]="SRR11674015,SRR11673972"
    ["OdNE"]="SRR11674008,SRR11673970"
    ["OdNW"]="SRR11674005,SRR11673969"
)

# CRISPR regions: Contig Start End ID
REGIONS=(
    "EsbE_1_2_3_4_6 368210 368412 CRISPR2341"
    "Hjor_1 38057 38746 CRISPR2344"
    "Hjor_1 81160 81847 CRISPR2345"
    "Hjor_1 355908 356129 CRISPR2346"
    "Hjor_1 377615 381005 CRISPR2347"
    "Lyne_1 263777 265136 CRISPR2348"
    "Lyne_1 274527 275508 CRISPR2349"
    "Lyne_1 422822 424179 CRISPR2350"
    "Lyne_1 628779 630646 CRISPR2351"
    "Lyne_1 1025331 1027131 CRISPR2352"
    "OdNE_1 18548 18756 CRISPR2353"
    "OdNE_1 785413 786062 CRISPR2354"
    "OdNE_1 861349 861600 CRISPR2355"
    "OdNE_1 1241114 1241881 CRISPR2356"
    "OdNE_1 1326127 1327053 CRISPR2357"
    "OdNE_1 1469951 1470166 CRISPR2358"
    "OdNW_1_7 64468 67523 CRISPR2359"
    "OdNW_1_7 86929 87856 CRISPR2360"
    "OdNW_1_7 88407 89639 CRISPR2361"
    "OdNW_1_7 148116 148977 CRISPR2362"
    "OdNW_1_7 171556 172476 CRISPR2363"
    "OdNW_1_7 273050 274662 CRISPR2364"
    "OdNW_1_7 317331 320793 CRISPR2365"
    "OdNW_1_7 339556 340227 CRISPR2366"
    "OdNW_1_7 399106 401326 CRISPR2367"
)

echo "========================================="
echo "Visualizing ${#REGIONS[@]} CRISPR regions"
echo "========================================="

# Counter
COUNT=0

# Loop through regions
for REGION in "${REGIONS[@]}"; do
    # Parse region info
    read -r CONTIG START END ID <<< "$REGION"
    
    # Extract station name (first part before underscore)
    STATION=$(echo "$CONTIG" | cut -d'_' -f1)
    
    # Get SRR accessions for this station
    SRRS="${SRR_MAP[$STATION]}"
    if [ -z "$SRRS" ]; then
        echo "ERROR: No SRR mapping found for station: $STATION"
        continue
    fi
    
    # Split into Illumina and ONT
    ILL_SRR=$(echo "$SRRS" | cut -d',' -f1)
    ONT_SRR=$(echo "$SRRS" | cut -d',' -f2)
    
    # Construct file paths
    FASTA="${FASTA_DIR}/${CONTIG}_polished.fasta"
    ONT_BAM="${ONT_DIR}/${ONT_SRR}_mapped_sorted.bam"
    ILL_BAM="${ILL_DIR}/${ILL_SRR}_mapped_sorted.bam"
    OUTPUT="${OUT_DIR}/${CONTIG}_${ID}_${START}-${END}.png"
    
    # Check if files exist
    if [ ! -f "$FASTA" ]; then
        echo "WARNING: FASTA not found: $FASTA - Skipping"
        continue
    fi
    if [ ! -f "$ONT_BAM" ]; then
        echo "WARNING: ONT BAM not found: $ONT_BAM - Skipping"
        continue
    fi
    if [ ! -f "$ILL_BAM" ]; then
        echo "WARNING: Illumina BAM not found: $ILL_BAM - Skipping"
        continue
    fi
    
    # Increment counter
    COUNT=$((COUNT + 1))
    
    # Run visualization
    echo ""
    echo "[$COUNT/${#REGIONS[@]}] Processing: ${CONTIG} ${START}-${END} (${ID})"
    
    Rscript scripts/61_visualize_array_region.R \
        "$FASTA" \
        "$ONT_BAM" \
        "$ILL_BAM" \
        "$CONTIG" \
        "$START" \
        "$END" \
        "$OUTPUT"
    
    echo "  -> Saved: $OUTPUT"
done

echo ""
echo "========================================="
echo "Complete! Processed $COUNT regions"
echo "Output directory: $OUT_DIR"
echo "========================================="
