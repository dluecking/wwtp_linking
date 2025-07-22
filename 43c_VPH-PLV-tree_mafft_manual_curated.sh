#!/bin/bash

module load mafft

mafft --auto intermediate/treebuilding/all_public_and_my_own_MCP_proteins_manually_cleaned.faa > intermediate/treebuilding/all_public_and_my_own_MCP_proteins_cleaned.aln
