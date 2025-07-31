#! /bin/bash

module load conda
conda activate aster

# concat the nuphylo trees
cat intermediate/GV_tree/nuhpylo_output/*/allseqs.nwk > intermediate/GV_tree/combined_nuphylo.nwk

# run the astral command
astral -i intermediate/GV_tree/combined_nuphylo.nwk -o intermediate/GV_tree/combined_nuphylo_astral.tre

