#!/bin/bash
#SBATCH --job-name=hmm_search
#SBATCH --output=log/output_hmm_search_%A_%a.out # OUTPUT THAT WOULD PRINT TO SCREEN GOES HERE
#SBATCH --error=log/error_hmm_search_%A_%a.err # LOG FILE GOES HERE
#SBATCH --array=1-23 # OPTIONAL! BUT OFTEN USEFUL
#SBATCH --ntasks=1 # FOR PARRALELIZATION, CAN SAFELY KEEP AT 1
#SBATCH --cpus-per-task=16 # NUMBER OF CPUS YOU NEED (good number: 4-16)
#SBATCH --time=5:00:00  # TIME YOU NEED DAY:HOURS:MINUTES:SECONDS (MORE -> it finishes probably, LESS -> you get priority on SLURM)
#SBATCH --mem=64G         # AMOUNT OF RAM YOU NEED (good number: 8GB for small jobs, 64GB for medium jobs, 256+GB for assemblies)

# load your modules here, e.g.:
module load hmmer
echo "loaded modules"

# get protein file
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" helperfiles/list_of_samples.txt)

PROTEIN_FILE=${line}.faa

hmmsearch \
  --cpu 16 \
  --tblout intermediate/hmm_out/${line}.tbl \
  --domtblout intermediate/hmm_out/${line}.domtbl \
  data/hmm/combined_vph_plv.hmm \
  data/proteins/${PROTEIN_FILE} \
  > intermediate/hmm_out/${line}.txt
