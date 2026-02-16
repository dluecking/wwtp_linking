#!/bin/bash
#SBATCH --job-name=34_dereplicate_lcs
#SBATCH --output=log/34_dereplicate_lcs_%j.log
#SBATCH --error=log/34_dereplicate_lcs_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G

# CONSTANTS and FILES
ALL_CONTIGS="intermediate/contigs/all_contigs.fna"
ALL_CONTIG_IDS="intermediate/network/all_contig_ids.txt"
PLV_VPH_GV_CONTIGS="intermediate/contigs/combined_gv_plv_vph.fna"
LC_CONTIGS="intermediate/contigs/lc_contigs_only.fna"
MASH_SKETCH="intermediate/mash/lc_contigs_only.msh"
MASH_DISTANCES="intermediate/mash/lc_contigs_only_mash_distances.txt"
MASH_DISTANCES_FILTERED="intermediate/mash/lc_contigs_only_mash_distances_filtered.txt"
CLUSTER_MEMBERSHIP="intermediate/mash/lc_contigs_only_cluster_membership.csv"

# threshold
CONTIG_DISSIMILARITY_THRESHOLD=0.05 # corresponds to 95% similarity

# load conda
module load Conda

### Part 0: only keep LC contigs for 
# load conda
conda activate bioinf
# Get non-LC IDs
grep "^>" "$PLV_VPH_GV_CONTIGS" | sed 's/^>//' > intermediate/mash/exclude_ids.txt
seqkit grep -n -v -f intermediate/mash/exclude_ids.txt "$ALL_CONTIGS" > "$LC_CONTIGS"

### PART 1 run mash sketch and mash distance
# load conda
conda activate mash-2.3

echo "Starting Mash sketching and distance calculation..."
mash sketch -i -p 16 -s 10000 -o ${MASH_SKETCH%.msh} $LC_CONTIGS
echo "Mash sketching complete. Now calculating pairwise distances..."
mash dist -p 16 $MASH_SKETCH $MASH_SKETCH > $MASH_DISTANCES
echo "Mash distance calculation complete. Now filtering distances based on threshold..."

### PART 2 simplify mash output
awk -v threshold=$CONTIG_DISSIMILARITY_THRESHOLD '$3 < threshold' $MASH_DISTANCES > $MASH_DISTANCES_FILTERED
echo "Filtering complete. Now running R script to cluster LCs based on Mash distances..."

### PART 3 run R script to dereplicate LCs based on mash output
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R_parallel
# Use ALL_CONTIG_IDS (not filtered to LC only)
Rscript scripts/34_cluster_lcs_based_on_mash_distances.R \
    --threshold $CONTIG_DISSIMILARITY_THRESHOLD \
    --mash_distances_filtered "$MASH_DISTANCES_FILTERED" \
    --all_lc_contig_ids "$ALL_CONTIG_IDS" \
    --membership_output "$CLUSTER_MEMBERSHIP"


echo "all done!"



