#!/bin/bash
#SBATCH --job-name=vcontact2_restart
#SBATCH --output=log/vcontact2_restart_%j.out
#SBATCH --error=log/vcontact2_restart_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    # Can try 16 first, then increase if needed
#SBATCH --time=3-00:00:00     # Start with 3 days, consider 5 days for buffer
#SBATCH --mem=512G             # Start with 64GB, consider 128GB if OOM occurs

module load conda
conda activate vcontact2

# make outdir
mkdir -p intermediate/vcontact2_out

# run the command
# vcontact2 --raw-proteins intermediate/proteins/all_proteins.faa \
#  --proteins-fp intermediate/proteins/all_proteins_gene_to_genome.csv \
#  --db 'None' --output-dir intermediate/vcontact2_out --c1-bin ~/bin/cluster_one-1.0.jar

# restarting the run
vcontact2 \
    --pcs intermediate/vcontact2_out/vConTACT_pcs.csv \
    --contigs intermediate/vcontact2_out/vConTACT_contigs.csv \
    --proteins-fp intermediate/proteins/all_proteins_gene_to_genome.csv \
    --pc-profiles intermediate/vcontact2_out/vConTACT_profiles.csv \
    --db 'None' \
    --output-dir intermediate/vcontact2_out \
    --c1-bin ~/bin/cluster_one-1.0.jar 


