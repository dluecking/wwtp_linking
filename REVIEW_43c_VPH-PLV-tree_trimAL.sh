#!/bin/bash

# manually curating is not accepted anymore
# doing the trimming of the original alignment with trimAl instead


module load trimAl/1.5.1-GCCcore-13.3.0

trimal \
  -in intermediate/treebuilding/all_public_and_my_own_MCP_proteins.aln \
  -out intermediate/treebuilding/all_public_and_my_own_MCP_proteins_trimal.aln \
  -automated1
  
