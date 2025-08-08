#!/bin/bash
# assign_taxonomy.sh
# Usage: ./assign_taxonomy.sh [--tmp <tmpdir>] [--threads <n>] input.fasta output.csv
#
# This script:
# 1. Uses a user-defined temporary directory (if provided) or creates one.
# 2. Runs Prodigal (--meta mode) to predict proteins.
# 3. Runs DIAMOND BLASTP on the predicted proteins against a Diamond-formatted NR database.
# 4. Parses DIAMOND output to extract, for each protein, the best hit's organism name (from brackets) and taxid.
# 5. Condenses per-protein hits to a per-contig taxonomic assignment, outputting:
#    contig_id, majority_scientific_name, majority_taxid, pct_majority_hits
#
# Required modules: prodigal, diamond

# Load required modules
module load prodigal
module load diamond

# Default settings
CUSTOM_TMPDIR=""
THREADS=16
POSITIONAL_ARGS=()

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --tmp)
            CUSTOM_TMPDIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Restore positional arguments
set -- "${POSITIONAL_ARGS[@]}"

# Check that exactly 2 positional arguments remain: input_fasta and output_csv.
if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 [--tmp <tmpdir>] [--threads <n>] input_fasta output_csv"
    exit 1
fi

INPUT_FASTA="$1"
OUTPUT_CSV="$2"

# Determine temporary working directory
if [[ -n "$CUSTOM_TMPDIR" ]]; then
    TMPDIR="$CUSTOM_TMPDIR"
    mkdir -p "$TMPDIR"
else
    TMPDIR=$(mktemp -d -t taxassign-XXXXXX)
fi

echo "Working in temporary directory: $TMPDIR"

# Step 0: Copy input.fasta if it doesn't exist in TMPDIR.
if [[ ! -f "$TMPDIR/input.fasta" ]]; then
    echo "Copying input FASTA to TMPDIR..."
    cp "$INPUT_FASTA" "$TMPDIR/input.fasta"
else
    echo "input.fasta already exists, skipping copy."
fi

# original directory
ORIG_DIR=$(pwd)

# move to tmpdir
cd "$TMPDIR" || exit 1

# Step 1: Run Prodigal (--meta mode) to predict proteins (output: proteins.faa).
if [[ ! -f "proteins.faa" ]]; then
    echo "Running Prodigal..."
    prodigal -i input.fasta -a proteins.faa -p meta
else
    echo "proteins.faa exists, skipping Prodigal."
fi

# copy nr_diamond to TMPDIR for speed up purposes
echo "Copying nr_diamond_tax.dmnd to $TMPDIR..."

cp /lisc/scratch/dome/luecking/dbs/diamond/nr_diamond_tax.dmnd "$TMPDIR/nr_diamond_tax.dmnd"
echo "Copied nr_diamond_tax.dmnd to $TMPDIR"

# Step 2: Run DIAMOND BLASTP against the NR database (output: diamond.out).
if [[ ! -f "diamond.out" ]]; then
    echo "Running DIAMOND BLASTP with $THREADS threads..."
    diamond blastp -q proteins.faa \
        -d "$TMPDIR/nr_diamond_tax.dmnd" \
        -o diamond.out \
        -f 6 qseqid stitle staxids \
        --threads "$THREADS" \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --faster \
        --tmpdir "$TMPDIR"
else
    echo "diamond.out exists, skipping DIAMOND BLASTP."
fi

# Step 3: Parse protein headers to count total predicted proteins per contig (output: total_proteins.txt).
if [[ ! -f "total_proteins.txt" ]]; then
    echo "Generating total_proteins.txt..."
    grep "^>" proteins.faa | sed 's/^>//' | awk '{
        match($1, /^(.*)_([0-9]+)$/, arr);
        if (arr[1] != "") {
            print arr[1];
        }
    }' | sort | uniq -c | awk '{print $2 "," $1}' > total_proteins.txt
else
    echo "total_proteins.txt exists, skipping this step."
fi

# Step 4: Parse DIAMOND output to extract contig, organism name (from brackets), and taxid (output: blast_parsed.txt).
if [[ ! -f "blast_parsed.txt" ]]; then
    echo "Parsing DIAMOND output to blast_parsed.txt..."
    awk '{
        # Remove last underscore and trailing digits from query ID
        contig = $1
        sub(/_[0-9]+$/, "", contig)

        # Rebuild stitle from fields 2 to end
        stitle = "";
        for (i = 2; i <= NF; i++) {
            stitle = stitle $i " ";
        }

        # Extract organism name from brackets [ ]
        match(stitle, /\[([^]]+)\]/, org);
        org_name = (org[1] != "") ? org[1] : "Unknown";

        print contig "\t" org_name;
    }' diamond.out > blast_parsed.txt
else
    echo "blast_parsed.txt exists, skipping parsing DIAMOND output."
fi

# Step 5: Compute per-contig taxonomy assignment (output: per_contig_taxonomy.csv).
if [[ ! -f "per_contig_taxonomy.csv" ]]; then
    echo "Generating summary per contig to per_contig_taxonomy.csv..."
    awk -F"," -v OFS="," '
        FNR==NR {
            total[$1] = $2;
            next;
        }

        # For second file, set FS to tab
        FNR==1 { FS="\t" }

        {
            contig = $1;
            organism = $2;
            blast_count[contig]++;
            hit[contig, organism]++;
        }

        END {
            print "contig_id","majority_organism","pct_reads_assigned_to_majority_taxon";
            for (c in total) {
                bc = blast_count[c] + 0;
                if (bc == 0) {
                    print c, "NA", 0;
                } else {
                    max = 0; majority = "";
                    for (k in hit) {
                        split(k, arr, SUBSEP);
                        if (arr[1] == c && hit[k] > max) {
                            max = hit[k];
                            majority = arr[2];
                        }
                    }
                    pct = (max / bc) * 100;
                    print c, majority, pct;
                }
            }
        }
    ' total_proteins.txt blast_parsed.txt > per_contig_taxonomy.csv
else
    echo "per_contig_taxonomy.csv exists, skipping summary step."
fi

# Move the result to the specified output CSV file (use original directory)
echo "Before moving: here is the current directory:"
ls -l
echo "Moving per_contig_taxonomy.csv to $ORIG_DIR/$OUTPUT_CSV"
cp per_contig_taxonomy.csv "$ORIG_DIR/$OUTPUT_CSV"
echo "Result saved to $ORIG_DIR/$OUTPUT_CSV"

# Note: We do not remove the temporary directory if a custom directory was provided.
if [[ -z "$CUSTOM_TMPDIR" ]]; then
    # Uncomment the following line to clean up auto-generated TMPDIR:
    # rm -rf "$TMPDIR"
    echo "Auto-generated temporary directory $TMPDIR remains (cleanup commented out)."
else
    echo "Custom TMPDIR provided; not removing $TMPDIR."
fi