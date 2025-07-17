#!/bin/bash


# STEP 1: Run Prodigal to predict genes of public vphs/plvs
echo "Running prodigal..."
module load prodigal

prodigal -i data/virophages_public/all_public_vphs.fasta -a intermediate/proteins/all_public_vphs_proteins.faa -p meta
echo "Prodigal completed."

# STEP 2: Run HMMER to search for MCPs
echo "Running hmmsearch..."
module load hmmer

PROTEIN_FILE=intermediate/proteins/all_public_vphs_proteins.faa
hmmsearch \
  --cpu 8 \
  -T 50 \
  --tblout intermediate/hmm_out/all_public_vphs.tbl \
  --domtblout intermediate/hmm_out/all_public_vphs.domtbl \
  data/hmm/combined_vph_plv.hmm \
  ${PROTEIN_FILE} \
  > intermediate/hmm_out/all_public_vphs.txt

echo "HMMER search completed."

# # STEP 3: Retrieve sequences of MCPs
# echo "Retrieving MCP sequences..."
# module load conda
# conda activate R_new

# Rscript scripts/retrieve_mcp_sequences.R  \
#   --public_hmmout intermediate/hmm_out/all_public_vphs.tbl \
#   --public_proteins intermediate/proteins/all_public_vphs_proteins.faa \
#   --my_hmmout_dir intermediate/hmm_out \
#   --my_vph_protein_dir intermediate/proteins/cleaned/vph \
#   --my_plv_protein_dir intermediate/proteins/cleaned/plv \
#   --output intermediate/proteins/all_public_and_my_own_MCP_proteins.faa 



# echo "MCP sequences retrieved."




# # STEP 4: Align sequences using MAFFT
# echo "Running MAFFT..."
# module load mafft

# mafft --auto intermediate/proteins/all_public_vphs_proteins_mcp.faa > intermediate/alignments/all_public_vphs_proteins_mcp_aligned.faa

# echo "MAFFT alignment completed."
