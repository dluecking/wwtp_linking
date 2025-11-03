import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import click # Import the click library

# CONSTANTS
MIN_READ_LENGTH = 1000  # Minimum total read length for a read to be considered
MIN_MAPPED_LENGTH = 500  # Minimum mapped length for a read to be considered
MIN_OVERHANG_LENGTH = 100  # Minimum length of soft-clipped overhang to be extracted
SOFTCLIP_WINDOW = 50  # Window size (in bp) from reference start/end to consider an alignment "at the boundary"

@click.command()
@click.option('--input-bam',
              type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True),
              required=True,
              help='Path to the input sorted BAM file.')
@click.option('--fasta-output-dir',
              type=click.Path(file_okay=False, writable=True, resolve_path=True),
              required=True,
              help='Directory to save output FASTA files. Will be created if it does not exist.')
@click.option('--bam-output-dir',
              type=click.Path(file_okay=False, writable=True, resolve_path=True),
              required=True,
              help='Directory to save output BAM files. Will be created if it does not exist.')
def extract_overhangs(input_bam, fasta_output_dir, bam_output_dir):
    """
    Extracts overhanging reads and their sequences from a BAM file,
    categorizing them by their alignment to reference contigs.
    Output BAM and FASTA files are generated per reference into specified directories.
    """
    click.echo(f"Starting overhang extraction with parameters:")
    click.echo(f"  Input BAM: {input_bam}")
    click.echo(f"  FASTA Output Directory: {fasta_output_dir}")
    click.echo(f"  BAM Output Directory: {bam_output_dir}")
    click.echo(f"  Min Read Length: {MIN_READ_LENGTH}")
    click.echo(f"  Min Mapped Length: {MIN_MAPPED_LENGTH}")
    click.echo(f"  Min Overhang Length: {MIN_OVERHANG_LENGTH}")
    click.echo(f"  Softclip Window: {SOFTCLIP_WINDOW}")

    # Ensure the output directories exist
    try:
        os.makedirs(fasta_output_dir, exist_ok=True)
        click.echo(f"Ensured FASTA output directory exists: {fasta_output_dir}")
        os.makedirs(bam_output_dir, exist_ok=True)
        click.echo(f"Ensured BAM output directory exists: {bam_output_dir}")
    except Exception as e:
        click.echo(f"Error creating output directories: {e}", err=True)
        return # Exit the command if directories cannot be created

    # Open BAM input file
    try:
        bamfile = pysam.AlignmentFile(input_bam, "rb")
        click.echo(f"Successfully opened input BAM file: {input_bam}")
    except FileNotFoundError:
        click.echo(f"Error: Input BAM file not found at {input_bam}. Please check the path.", err=True)
        return
    except Exception as e:
        click.echo(f"Error opening input BAM file {input_bam}: {e}", err=True)
        return

    # Dictionaries to store output file handles and FASTA records, keyed by reference name
    output_bam_files = {}
    fasta_out_per_ref = {}

    # Store reference lengths in a dictionary for quick lookup
    ref_lengths = dict(zip(bamfile.references, bamfile.lengths))
    click.echo(f"Found {len(ref_lengths)} references in BAM header.")

    total_processed_reads_count = 0

    # Iterate over each read in the BAM file
    # and handles unsorted BAMs without requiring an index for full iteration.
    for read in bamfile.fetch(until_eof=True):
        read_id = read.query_name
        # Using click.echo directly to avoid issues with standard print in some environments
        click.echo(f"\n--- Processing read: {read_id} ---")

        # Skip unmapped reads
        if read.is_unmapped:
            click.echo(f"  -> Discarded: Read {read_id} is unmapped.")
            continue

        # Get mapped length and total read length
        mapped_len = read.query_alignment_length
        total_read_len = read.query_length
        click.echo(f"  -> Total read length: {total_read_len} bp, Mapped length: {mapped_len} bp")

        # Skip reads that don't meet minimum length criteria
        if total_read_len is None or mapped_len < MIN_MAPPED_LENGTH or total_read_len < MIN_READ_LENGTH:
            click.echo(f"  -> Discarded: Read {read_id} fails length criteria (Min total: {MIN_READ_LENGTH}, Min mapped: {MIN_MAPPED_LENGTH}).")
            continue

        # Get reference name, reference length, and mapped start/end positions
        rname = read.reference_name
        ref_len = ref_lengths.get(rname)
        if ref_len is None:
            click.echo(f"Warning: Reference '{rname}' not found in header for read {read_id}. Skipping.", err=True)
            continue

        start = read.reference_start  # 0-based start position on reference
        end = read.reference_end      # 0-based end position (exclusive) on reference
        click.echo(f"  -> Mapped to '{rname}' from {start} to {end} (Reference length: {ref_len} bp)")

        # Get CIGAR tuples. pysam returns a list of (operation, length) tuples.
        cigars = read.cigartuples or []
        click.echo(f"  -> CIGAR tuples: {cigars}")

        # --- Determine actual soft-clip lengths from CIGAR at the ends of the read ---
        overhang_5_actual_len = 0
        if cigars and cigars[0][0] == 4: # Check if the first CIGAR operation is Soft Clip (S)
            overhang_5_actual_len = cigars[0][1]
        click.echo(f"  -> Raw 5' soft clip length: {overhang_5_actual_len} bp")

        overhang_3_actual_len = 0
        if cigars and cigars[-1][0] == 4: # Check if the last CIGAR operation is Soft Clip (S)
            overhang_3_actual_len = cigars[-1][1]
        click.echo(f"  -> Raw 3' soft clip length: {overhang_3_actual_len} bp")

        # --- Classify potential overhangs based on the new definitions ---

        # Condition for a 5' soft clip that *extends out* of the reference boundary
        is_5prime_extending_overhang = (
            overhang_5_actual_len >= MIN_OVERHANG_LENGTH and
            start <= SOFTCLIP_WINDOW
        )
        click.echo(f"  -> Qualifies for _overhang_5 (5' extending): (Len >= {MIN_OVERHANG_LENGTH}? {overhang_5_actual_len >= MIN_OVERHANG_LENGTH}) AND (Start <= {SOFTCLIP_WINDOW}? {start <= SOFTCLIP_WINDOW}) -> {is_5prime_extending_overhang}")

        # Condition for a 3' soft clip that *extends out* of the reference boundary
        is_3prime_extending_overhang = (
            overhang_3_actual_len >= MIN_OVERHANG_LENGTH and
            end >= (ref_len - SOFTCLIP_WINDOW)
        )
        click.echo(f"  -> Qualifies for _overhang_3 (3' extending): (Len >= {MIN_OVERHANG_LENGTH}? {overhang_3_actual_len >= MIN_OVERHANG_LENGTH}) AND (End >= (RefLen-{SOFTCLIP_WINDOW})? {end >= (ref_len - SOFTCLIP_WINDOW)}) -> {is_3prime_extending_overhang}")

        # Condition for a 5' soft clip that is *internal* to the reference alignment
        is_5prime_internal_softclip = (
            overhang_5_actual_len >= MIN_OVERHANG_LENGTH and
            start > SOFTCLIP_WINDOW
        )
        click.echo(f"  -> Qualifies for _overhang_m (5' internal soft clip): (Len >= {MIN_OVERHANG_LENGTH}? {overhang_5_actual_len >= MIN_OVERHANG_LENGTH}) AND (Start > {SOFTCLIP_WINDOW}? {start > SOFTCLIP_WINDOW}) -> {is_5prime_internal_softclip}")

        # Condition for a 3' soft clip that is *internal* to the reference alignment
        is_3prime_internal_softclip = (
            overhang_3_actual_len >= MIN_OVERHANG_LENGTH and
            end < (ref_len - SOFTCLIP_WINDOW)
        )
        click.echo(f"  -> Qualifies for _overhang_m (3' internal soft clip): (Len >= {MIN_OVERHANG_LENGTH}? {overhang_3_actual_len >= MIN_OVERHANG_LENGTH}) AND (End < (RefLen-{SOFTCLIP_WINDOW})? {end < (ref_len - SOFTCLIP_WINDOW)}) -> {is_3prime_internal_softclip}")

        # Decide if the read should be written to BAM based on *any* qualifying overhang type
        should_write_to_bam = (
            is_5prime_extending_overhang or
            is_3prime_extending_overhang or
            is_5prime_internal_softclip or
            is_3prime_internal_softclip
        )

        if not should_write_to_bam:
            click.echo(f"  -> Discarded: Read {read_id} has no qualifying overhangs (min length {MIN_OVERHANG_LENGTH}bp) based on current definitions.")
            continue

        # --- Write the read to the appropriate output BAM file ---
        if rname not in output_bam_files:
            output_bam_filename = os.path.join(bam_output_dir, f"overhanging_{rname}.bam")
            try:
                output_bam_files[rname] = pysam.AlignmentFile(output_bam_filename, "wb", template=bamfile)
                click.echo(f"  -> Created new BAM output file: {output_bam_filename}")
            except Exception as e:
                click.echo(f"Error creating BAM output file {output_bam_filename}: {e}", err=True)
                # If we can't create the BAM, we should probably skip processing FASTA for this read too
                continue

        output_bam_files[rname].write(read)
        click.echo(f"  -> Read {read_id} written to output BAM file: {os.path.basename(output_bam_filename)}")
        total_processed_reads_count += 1

        # --- Prepare for FASTA extraction ---
        seq = read.query_sequence
        is_rev = read.is_reverse # True if the read is mapped to the reverse strand
        click.echo(f"  -> Read is_reverse: {is_rev}")

        # Initialize FASTA list for this reference if not already done
        if rname not in fasta_out_per_ref:
            fasta_out_per_ref[rname] = []

        # --- Extract and write FASTA entries based on their specific classifications ---
        # Each condition is checked independently to capture all applicable overhangs for a read.

        if is_5prime_extending_overhang:
            over_seq = seq[:overhang_5_actual_len]
            if is_rev:
                over_seq = str(Seq(over_seq).reverse_complement())
            fasta_out_per_ref[rname].append(SeqRecord(Seq(over_seq), id=f"{read_id}_to_{rname}_overhang_5", description=f"5' boundary overhang on {rname}"))
            click.echo(f"    -> Added FASTA entry: {read_id}_to_{rname}_overhang_5 (5' boundary, length {overhang_5_actual_len}bp) for {rname}")

        if is_3prime_extending_overhang:
            over_seq = seq[-overhang_3_actual_len:]
            if is_rev:
                over_seq = str(Seq(over_seq).reverse_complement())
            fasta_out_per_ref[rname].append(SeqRecord(Seq(over_seq), id=f"{read_id}_to_{rname}_overhang_3", description=f"3' boundary overhang on {rname}"))
            click.echo(f"    -> Added FASTA entry: {read_id}_to_{rname}_overhang_3 (3' boundary, length {overhang_3_actual_len}bp) for {rname}")

        if is_5prime_internal_softclip:
            over_seq = seq[:overhang_5_actual_len]
            if is_rev:
                over_seq = str(Seq(over_seq).reverse_complement())
            fasta_out_per_ref[rname].append(SeqRecord(Seq(over_seq), id=f"{read_id}_to_{rname}_overhang_m", description=f"5' internal soft clip on {rname}"))
            click.echo(f"    -> Added FASTA entry: {read_id}_to_{rname}_overhang_m (5' internal, length {overhang_5_actual_len}bp) for {rname}")

        if is_3prime_internal_softclip:
            over_seq = seq[-overhang_3_actual_len:]
            if is_rev:
                over_seq = str(Seq(over_seq).reverse_complement())
            fasta_out_per_ref[rname].append(SeqRecord(Seq(over_seq), id=f"{read_id}_to_{rname}_overhang_m", description=f"3' internal soft clip on {rname}"))
            click.echo(f"    -> Added FASTA entry: {read_id}_to_{rname}_overhang_m (3' internal, length {overhang_3_actual_len}bp) for {rname}")


    # --- Write all collected overhang sequences to their respective FASTA files ---
    click.echo("\n--- Writing FASTA output files ---")
    for rname, fasta_records in fasta_out_per_ref.items():
        if fasta_records: # Only write if there are records for this reference
            fasta_filename = os.path.join(fasta_output_dir, f"overhangs_{rname}.fasta")
            try:
                SeqIO.write(fasta_records, fasta_filename, "fasta")
                click.echo(f"Successfully wrote {len(fasta_records)} overhang sequences to {fasta_filename}")
            except IOError as e:
                click.echo(f"Error writing FASTA file {fasta_filename}: {e}", err=True)
                click.echo("FASTA output might be empty or problematic.", err=True)
        else:
            click.echo(f"No overhangs found for reference '{rname}', skipping FASTA file creation.")


    # --- Close all opened BAM files ---
    click.echo("\n--- Closing BAM output files ---")
    for rname, outfile_handle in output_bam_files.items():
        try:
            outfile_handle.close()
            click.echo(f"Closed BAM output file for {rname}.")
        except Exception as e:
            click.echo(f"Error closing BAM file for {rname}: {e}", err=True)


    # Close the input BAM file
    bamfile.close()
    click.echo(f"\nProcessing complete. {total_processed_reads_count} reads were written to BAM files across all references. All files closed.")

if __name__ == '__main__':
    extract_overhangs()
