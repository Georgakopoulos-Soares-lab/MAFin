import sys
import re
import argparse
import json
import logging
import os
import multiprocessing
from Bio import AlignIO
from io import StringIO
from Bio import motifs
from Bio.Seq import Seq
import random

# Configure logging to write to 'logging.txt'
logging.basicConfig(
    filename='logging.txt',
    filemode='w',
    level=logging.INFO,
    format='%(message)s'
)

def load_genome_ids(genome_ids_file):
    """Load genome IDs from a given file.

    The genome IDs are used to specify which genomes to consider for conservation calculations.
    If not provided, all genomes in the MAF file will be considered.
    """
    if genome_ids_file:
        logging.info(f"Loading genome IDs from {genome_ids_file}")
        with open(genome_ids_file, 'r') as file:
            return file.read().splitlines()
    return None  # If not provided, we consider all genomes

def assign_file_chunks(maf_file, num_processes):
    """Assign chunks of the MAF file to processes based on byte positions."""
    file_size = os.path.getsize(maf_file)
    chunk_size = file_size // num_processes
    chunk_positions = []

    with open(maf_file, 'rb') as f:
        for i in range(num_processes):
            start_pos = i * chunk_size
            if start_pos > 0:
                f.seek(start_pos)
                # Move to the end of the current line
                f.readline()
            # Now, find the next line that starts with 'a'
            line = f.readline()
            while line and not line.startswith(b'a'):
                line = f.readline()
            # Record the position
            chunk_start = f.tell() - len(line)
            chunk_positions.append(chunk_start)

    # Now, for end positions, it's the start of the next chunk, or EOF for the last one
    chunk_ranges = []
    for i in range(num_processes):
        start_pos = chunk_positions[i]
        end_pos = chunk_positions[i + 1] if i + 1 < len(chunk_positions) else os.path.getsize(maf_file)
        chunk_ranges.append((start_pos, end_pos))

    return chunk_ranges

def parse_block_from_string(block_data):
    """Parse a MAF block from a string into a MultipleSeqAlignment object."""
    handle = StringIO(block_data)
    try:
        alignment = AlignIO.read(handle, "maf")
        return alignment
    except Exception as e:
        logging.error(f"Error parsing block: {e}", exc_info=True)
        return None

def get_reverse_complement_gapped_sequence(seq_gapped):
    """Reverse complement a gapped sequence."""
    complement_map = str.maketrans('ACGTacgt', 'TGCAtgca')
    reversed_seq = seq_gapped[::-1]
    reverse_complemented_seq = ''.join(
        base.translate(complement_map) if base != '-' else '-'
        for base in reversed_seq
    )
    return reverse_complemented_seq

def compute_threshold_score(pssm, pvalue, num_samples=10000):
    """Compute the threshold score for a given p-value by sampling random sequences."""
    logging.info(f"Computing threshold score for p-value: {pvalue}")
    scores = []
    for _ in range(num_samples):
        seq = generate_random_sequence(pssm.length, pssm.background)
        score = pssm.calculate(seq)
        scores.append(score)
    scores.sort()
    index = int((1 - pvalue) * num_samples)
    threshold_score = scores[index] if index < len(scores) else scores[-1]
    logging.info(f"Threshold score corresponding to p-value {pvalue}: {threshold_score}")
    return threshold_score

def generate_random_sequence(length, background):
    """Generate a random sequence based on the background nucleotide frequencies."""
    nucleotides = ['A', 'C', 'G', 'T']
    frequencies = [background.get(nuc, 0.25) for nuc in nucleotides]
    return Seq(''.join(random.choices(nucleotides, frequencies, k=length)))

def compute_vectors_and_conservation(block, genome_ids, gapped_start, gapped_end, ref_gapped_seq, ref_genome_id, ref_seq_record):
    """
    Compute similarity vectors and conservation percentages for the given gapped positions,
    based on the gapped motif sequence.
    """
    vectors = {}
    aligned_sequences = []
    total_conservation = 0
    genomes_with_data = 0
    motif_length = len(ref_gapped_seq.replace('-', ''))  # True length of the motif (excluding positions where both have gaps)

    # Compute ungapped positions before motif start in reference genome
    ref_seq_gapped = str(ref_seq_record.seq)
    ref_strand = ref_seq_record.annotations.get('strand', '+')
    ref_start = int(ref_seq_record.annotations.get('start', 0))
    ref_size = int(ref_seq_record.annotations.get('size', 0))
    ref_src_size = int(ref_seq_record.annotations.get('srcSize', 0))

    # Compute ungapped positions before motif start
    ungapped_positions_before_motif = len([c for c in ref_seq_gapped[:gapped_start] if c != '-'])
    motif_ungapped_length = len(ref_gapped_seq.replace('-', ''))

    # Compute genomic start and end positions for reference genome
    if ref_strand == '+':
        genomic_start = ref_start + ungapped_positions_before_motif
        genomic_end = genomic_start + motif_ungapped_length - 1
    else:
        genomic_end = ref_start + ref_size - ungapped_positions_before_motif
        genomic_start = genomic_end - motif_ungapped_length + 1

    # Extract chromosome information
    ref_chromosome = '.'.join(ref_seq_record.id.split('.')[1:])  # Skip the genome ID

    for seq_record in block:
        genome_id = seq_record.id.split('.')[0]
        if genome_ids is None or genome_id in genome_ids:
            seq_gapped = str(seq_record.seq)

            # Extract the corresponding gapped sequence fragment
            seq_gapped_fragment = seq_gapped[gapped_start:gapped_end]

            # Initialize vector and matches
            vector = ''
            matches = 0
            positions_processed = 0  # Number of positions processed (excluding positions where both have gaps)

            ref_seq_fragment = ref_gapped_seq
            seq_seq_fragment = seq_gapped_fragment

            i = 0  # Position index
            while positions_processed < motif_length and i < len(ref_seq_fragment):
                ref_base = ref_seq_fragment[i]
                seq_base = seq_seq_fragment[i] if i < len(seq_seq_fragment) else '-'

                if ref_base == '-' and seq_base == '-':
                    # Both have gaps; skip this position
                    i += 1
                    continue
                else:
                    # Compare the bases
                    if ref_base == seq_base:
                        vector += '1'
                        matches += 1
                    else:
                        vector += '0'
                    positions_processed += 1
                    i += 1

            # Ensure the vector length matches the motif length
            if len(vector) < motif_length:
                zeros_to_add = motif_length - len(vector)
                vector += '0' * zeros_to_add

            # Calculate conservation percentage
            conservation_pct = (matches / motif_length) * 100 if motif_length > 0 else 0.0

            if genome_id != ref_genome_id:
                total_conservation += conservation_pct
                genomes_with_data += 1

            # Compute genomic coordinates for this genome
            seq_strand = seq_record.annotations.get('strand', '+')
            seq_start = int(seq_record.annotations.get('start', 0))
            seq_size = int(seq_record.annotations.get('size', 0))
            seq_src_size = int(seq_record.annotations.get('srcSize', 0))
            seq_id_parts = seq_record.id.split('.')
            seq_chromosome = '.'.join(seq_id_parts[1:])  # Skip the genome ID

            # Compute ungapped positions before motif start in this genome
            seq_gapped_full = str(seq_record.seq)
            ungapped_positions_before_motif_seq = len([c for c in seq_gapped_full[:gapped_start] if c != '-'])
            motif_ungapped_length_seq = len(seq_gapped_fragment.replace('-', ''))

            if seq_strand == '+':
                seq_genomic_start = seq_start + ungapped_positions_before_motif_seq
                seq_genomic_end = seq_genomic_start + motif_ungapped_length_seq - 1
            else:
                seq_genomic_end = seq_start + seq_size - ungapped_positions_before_motif_seq
                seq_genomic_start = seq_genomic_end - motif_ungapped_length_seq + 1

            aligned_sequences.append({
                "genome_id": genome_id,
                "chromosome": seq_chromosome,
                "genomic_start": seq_genomic_start,
                "genomic_end": seq_genomic_end,
                "strand": seq_strand,
                "aligned_sequence_gapped": seq_gapped_fragment,
                "vector": vector,
                "conservation": f"{conservation_pct:.2f}%",
                "gapped_start": int(gapped_start + 1),
                "gapped_end": int(gapped_end),
            })

            vectors[genome_id] = vector

    # Calculate average conservation value
    conservation_value = (total_conservation / genomes_with_data) if genomes_with_data > 0 else 0.0

    return vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome

def search_patterns_in_block(block, search_type, patterns, genome_ids, block_no, search_in, genome_files, reverse_complement, process_id):
    """Search for motifs in the sequences based on the specified search type."""
    try:
        # Determine the reference genome in this block
        ref_seq_record = None
        if search_in == 'reference':
            # The reference genome is the first sequence in the block
            ref_seq_record = block[0]
            sequences_to_search = [ref_seq_record]
        elif search_in == 'all':
            # Will search in all sequences
            if genome_ids is None:
                sequences_to_search = block
            else:
                sequences_to_search = [seq for seq in block if seq.id.split('.')[0] in genome_ids]
        else:
            logging.error(f"Invalid value for search_in: {search_in}")
            sys.exit(1)

        for seq_record in sequences_to_search:
            ref_genome_id = seq_record.id.split('.')[0]

            # Create output files for genomes if not already created
            if ref_genome_id not in genome_files:
                genome_dir = os.path.join('results', ref_genome_id)
                os.makedirs(genome_dir, exist_ok=True)
                output_file = os.path.join(genome_dir, f'motif_hits_process_{process_id}.json')
                genome_file = open(output_file, 'a')
                genome_files[ref_genome_id] = genome_file

            genome_file = genome_files[ref_genome_id]  # Get the file handle for the genome

            # Extract sequence information
            seq_gapped = str(seq_record.seq)  # Gapped sequence
            sequence_length = len(seq_gapped)
            sequence_strand = seq_record.annotations.get('strand', '+')
            if sequence_strand == -1:
                sequence_strand = '-'
            elif sequence_strand == 1:
                sequence_strand = '+'
            else:
                sequence_strand = '+'

            # Prepare sequences for matching
            seq_gapped_orig = seq_gapped
            seq_ungapped_orig = seq_gapped_orig.replace('-', '')

            if search_type == 'pwm':
                pwm_motifs = patterns
                seq_obj = Seq(seq_ungapped_orig)

                for pwm_data in pwm_motifs:
                    pwm = pwm_data['pwm']
                    pwm_name = pwm_data['name']
                    threshold_score = pwm_data['threshold']

                    # Search in original sequence
                    for position in range(len(seq_obj) - pwm.length + 1):
                        window_seq = seq_obj[position:position + pwm.length]
                        score = pwm.calculate(window_seq)
                        if score >= threshold_score:
                            pos = position
                            motif_len = pwm.length
                            ungapped_sequence = str(window_seq)

                            # Map ungapped positions to gapped positions
                            ungapped_to_gapped = [i for i, c in enumerate(seq_gapped_orig) if c != '-']
                            gapped_start = ungapped_to_gapped[pos]
                            gapped_end = ungapped_to_gapped[pos + motif_len - 1] + 1

                            # Extract gapped sequence
                            gapped_sequence = seq_gapped_orig[gapped_start:gapped_end]

                            # Determine motif strand
                            motif_strand = sequence_strand

                            # Compute conservation and genomic coordinates
                            ref_gapped_seq = gapped_sequence
                            vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                                block, genome_ids, gapped_start, gapped_end, ref_gapped_seq, ref_genome_id, seq_record
                            )

                            # Collect motif hit data
                            motif_hit_data = {
                                "hit_id": f"pwm_{ref_genome_id}_{pwm_name}_{block_no}_{pos}",
                                "block_no": block_no,
                                "reference_genome": ref_genome_id,
                                "chromosome": ref_chromosome,
                                "genomic_start": genomic_start,
                                "genomic_end": genomic_end,
                                "strand": sequence_strand,
                                "motif_type": "pwm",
                                "motif_length": motif_len,
                                "motif": ungapped_sequence,
                                "gapped_start": gapped_start + 1,
                                "gapped_end": gapped_end,
                                "ungapped_start": pos + 1,
                                "ungapped_end": pos + motif_len,
                                "motif_strand": motif_strand,
                                "score": float(score),  # Convert to float
                                "conservation": {
                                    "value": float(conservation_value),  # Convert to float
                                    "other_genomes": aligned_sequences
                                }
                            }

                            # Write the motif hit data to the JSON file
                            write_motif_hit(genome_file, motif_hit_data)

                    # Reverse complement PWM and search if reverse_complement is 'yes'
                    if reverse_complement == 'yes':
                        pwm_rc = pwm.reverse_complement()
                        pwm_name_rc = pwm_name + '_rc' if not pwm_name.endswith('_rc') else pwm_name
                        for position in range(len(seq_obj) - pwm_rc.length + 1):
                            window_seq = seq_obj[position:position + pwm_rc.length]
                            score = pwm_rc.calculate(window_seq)
                            if score >= threshold_score:
                                pos = position
                                motif_len = pwm_rc.length
                                ungapped_sequence = str(window_seq)

                                # Map ungapped positions to gapped positions
                                ungapped_to_gapped = [i for i, c in enumerate(seq_gapped_orig) if c != '-']
                                gapped_start = ungapped_to_gapped[pos]
                                gapped_end = ungapped_to_gapped[pos + motif_len - 1] + 1

                                # Extract gapped sequence
                                gapped_sequence = seq_gapped_orig[gapped_start:gapped_end]

                                # Determine motif strand
                                motif_strand = '-' if sequence_strand == '+' else '+'

                                # Compute conservation and genomic coordinates
                                ref_gapped_seq = gapped_sequence
                                vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                                    block, genome_ids, gapped_start, gapped_end, ref_gapped_seq, ref_genome_id, seq_record
                                )

                                # Collect motif hit data
                                motif_hit_data = {
                                    "hit_id": f"pwm_{ref_genome_id}_{pwm_name_rc}_{block_no}_{pos}",
                                    "block_no": block_no,
                                    "reference_genome": ref_genome_id,
                                    "chromosome": ref_chromosome,
                                    "genomic_start": genomic_start,
                                    "genomic_end": genomic_end,
                                    "strand": sequence_strand,
                                    "motif_type": "pwm",
                                    "motif_length": motif_len,
                                    "motif": ungapped_sequence,
                                    "gapped_start": gapped_start + 1,
                                    "gapped_end": gapped_end,
                                    "ungapped_start": pos + 1,
                                    "ungapped_end": pos + motif_len,
                                    "motif_strand": motif_strand,
                                    "score": float(score),  # Convert to float
                                    "conservation": {
                                        "value": float(conservation_value),  # Convert to float
                                        "other_genomes": aligned_sequences
                                    }
                                }

                                # Write the motif hit data to the JSON file
                                write_motif_hit(genome_file, motif_hit_data)

    except Exception as e:
        logging.error(f"Error processing block {block_no}: {e}", exc_info=True)
        sys.exit(1)  # Stop the script

def write_motif_hit(genome_file, motif_hit_data):
    """Write a motif hit to the genome's JSON file."""
    # Write the JSON object
    json.dump(motif_hit_data, genome_file)
    genome_file.write('\n')  # Write newline to separate JSON objects

def process_file_chunk(maf_file, start_pos, end_pos, search_type, patterns, genome_ids, search_in, reverse_complement, process_id):
    genome_files = {}

    try:
        with open(maf_file, 'rb') as handle:
            # Seek to the starting position
            handle.seek(start_pos)
            current_position = handle.tell()
            line = handle.readline()
            block_no = 0  # Initialize block number
            while current_position < end_pos and line:
                # Read until we find a line starting with 'a'
                while line and not line.startswith(b'a'):
                    current_position = handle.tell()
                    line = handle.readline()
                    if current_position >= end_pos:
                        break
                if current_position >= end_pos or not line:
                    break
                # Start of a block
                block_data = line.decode('utf-8')
                current_position = handle.tell()
                line = handle.readline()
                while line and not line.startswith(b'a'):
                    block_data += line.decode('utf-8')
                    current_position = handle.tell()
                    line = handle.readline()
                    if current_position >= end_pos:
                        break
                # Parse the block
                block = parse_block_from_string(block_data)
                if block:
                    # Process the block
                    search_patterns_in_block(block, search_type, patterns, genome_ids, block_no, search_in, genome_files, reverse_complement, process_id)
                block_no += 1
    except Exception as e:
        logging.error(f"Error in process {process_id}: {e}", exc_info=True)
        sys.exit(1)
    finally:
        # Close all genome files
        for genome_file in genome_files.values():
            genome_file.close()

def merge_results(results_dir):
    """Merge results from all processes into per-genome JSON files."""
    # Get list of genome directories
    genome_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    for genome_id in genome_dirs:
        merged_data = []
        genome_dir = os.path.join(results_dir, genome_id)
        # List all process-specific files
        process_files = [f for f in os.listdir(genome_dir) if f.startswith('motif_hits_process_') and f.endswith('.json')]
        # Sort the files to maintain order
        process_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))  # Extract process ID for sorting
        for process_file in process_files:
            genome_file_path = os.path.join(genome_dir, process_file)
            if os.path.exists(genome_file_path):
                with open(genome_file_path, 'r') as f:
                    for line in f:
                        try:
                            data = json.loads(line)
                            merged_data.append(data)
                        except json.JSONDecodeError:
                            logging.error(f"Error decoding JSON from {genome_file_path}")
        if merged_data:
            # Write merged data to final result file
            with open(os.path.join(genome_dir, f'motif_hits.json'), 'w') as f:
                json.dump(merged_data, f, indent=4)
            # Optionally, remove the process-specific files
            for process_file in process_files:
                os.remove(os.path.join(genome_dir, process_file))
        else:
            logging.info(f"No data found for genome {genome_id}")

def main(args):
    logging.info("Starting...")

    maf_file = args.maf_file
    genome_ids_file = args.genome_ids
    search_in = args.search_in
    reverse_complement = args.reverse_complement
    pvalue_threshold = args.pvalue_threshold
    num_processes = args.processes

    # Create the /results directory if it doesn't exist
    results_dir = os.path.join(os.getcwd(), 'results')
    os.makedirs(results_dir, exist_ok=True)

    # Load the genome IDs to consider from genome_ids.txt, if provided
    genome_ids = load_genome_ids(genome_ids_file)
    # If genome_ids is None, we consider all genomes

    # Determine search type and load patterns
    search_type = None
    patterns = None

    if args.regexes:
        search_type = 'regex'
        patterns = args.regexes
    elif args.kmers:
        search_type = 'kmer'
        kmer_patterns = []
        with open(args.kmers, 'r') as f:
            kmers = [line.strip() for line in f if line.strip()]
        kmer_patterns = [(kmer, '+') for kmer in kmers]
        # Generate reverse complements if reverse_complement is 'yes'
        if reverse_complement == 'yes':
            rc_kmers = [(str(Seq(kmer).reverse_complement()), '-') for kmer in kmers]
            kmer_patterns.extend(rc_kmers)
        # Remove duplicates
        kmer_patterns = list(set(kmer_patterns))
        patterns = kmer_patterns
    elif args.jaspar_file:
        search_type = 'pwm'
        pwm_motifs = []
        with open(args.jaspar_file, 'r') as f:
            content = f.read()
            if not content.startswith('>'):
                logging.error("PWM file format not recognized. Only JASPAR format is accepted.")
                sys.exit(1)
            f.seek(0)
            for motif in motifs.parse(f, 'jaspar'):
                # Generate PSSM
                pwm = motif.counts.normalize(pseudocounts=0.1)
                pssm = pwm.log_odds()
                pssm.background = motif.background if motif.background else {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
                # Compute threshold score
                threshold_score = compute_threshold_score(pssm, pvalue_threshold)
                pwm_motifs.append({'pwm': pssm, 'name': motif.name, 'threshold': threshold_score})
                # Reverse complement PWM is handled in the search function
        patterns = pwm_motifs
    else:
        logging.error("Error: One of --regexes, --kmers, or --jaspar_file must be provided.")
        sys.exit(1)

    # Assign chunks of the file to processes
    chunk_ranges = assign_file_chunks(maf_file, num_processes)

    # Initialize multiprocessing
    processes = []
    for i, (start_pos, end_pos) in enumerate(chunk_ranges):
        p = multiprocessing.Process(target=process_file_chunk, args=(
            maf_file, start_pos, end_pos, search_type, patterns, genome_ids, search_in, reverse_complement, i))
        processes.append(p)
        p.start()

    # Wait for all processes to finish
    for p in processes:
        p.join()

    # Merge results from all processes
    merge_results(results_dir)

    logging.info("Completed processing.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Search uncompressed MAF files for specific motifs using regex patterns, K-mers, or PWMs.")

    parser.add_argument('--maf_file', required=True, help="Path to the uncompressed MAF file.")
    parser.add_argument('--genome_ids', help="Path to the genome IDs file (optional). If not provided, all genomes are used.")
    parser.add_argument('--search_in', choices=['reference', 'all'], default='reference',
                        help="Specify whether to search motifs in 'reference' genome or 'all' genomes. Default is 'reference'.")
    parser.add_argument('--reverse_complement', choices=['yes', 'no'], default='no',
                        help="Specify whether to search motifs on both strands. Default is 'no'.")
    parser.add_argument('--pvalue_threshold', type=float, default=1e-4,
                        help="P-value threshold for PWM matches. Default is 1e-4.")
    parser.add_argument('--processes', type=int, default=1,
                        help="Number of parallel processes to use. Default is 1.")

    # Mutually exclusive group for search types
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--regexes', nargs='+', help="Regex patterns to search for in the MAF file.")
    group.add_argument('--kmers', help="Path to the K-mers file (one per line).")
    group.add_argument('--jaspar_file', help="Path to the PWM file in JASPAR format.")

    args = parser.parse_args()

    main(args)

