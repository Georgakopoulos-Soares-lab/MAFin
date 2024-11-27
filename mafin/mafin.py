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
import numpy as np
import csv
import shutil
import ahocorasick  # Updated import statement

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Configure logging to write to 'logging.txt'
logging.basicConfig(
    level=logging.INFO,  # Set logging level
    format='%(asctime)s - %(levelname)s - %(message)s',  # Set the format for logs
    handlers=[
        logging.FileHandler('logging.txt', mode='w'),  # Log to file
        logging.StreamHandler(sys.stdout)  # Log to stdout
    ]
)


def load_genome_ids(genome_ids_file):
    """Load genome IDs from a given file.

    The genome IDs are used to specify which genomes to consider for conservation calculations.
    If not provided, all genomes in the MAF file will be considered.
    """
    if genome_ids_file:
        logging.info(f"Loading genome IDs from {genome_ids_file}")
        try:
            with open(genome_ids_file, 'r') as file:
                genome_ids = [line.strip() for line in file if line.strip()]
            logging.info(f"Loaded {len(genome_ids)} genome IDs.")
            return genome_ids
        except Exception as e:
            logging.error(f"Error loading genome IDs from {genome_ids_file}: {e}", exc_info=True)
            sys.exit(1)
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
                f.readline()  # Move to the end of the current line
            chunk_start = f.tell()
            chunk_positions.append(chunk_start)

    # Define end positions
    chunk_ranges = []
    for i in range(num_processes):
        start_pos = chunk_positions[i]
        end_pos = chunk_positions[i + 1] if i + 1 < len(chunk_positions) else file_size
        chunk_ranges.append((start_pos, end_pos))

    logging.info(f"Assigned file chunks: {chunk_ranges}")
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
    logging.info(f"Computing threshold score for p-value: {pvalue} with background frequencies: {pssm.background}")
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
    total_freq = sum(frequencies)
    if not np.isclose(total_freq, 1.0):
        # Normalize frequencies
        frequencies = [freq / total_freq for freq in frequencies]
    return Seq(''.join(random.choices(nucleotides, frequencies, k=length)))

def compute_vectors_and_conservation(block, genome_ids, gapped_start, gapped_end,
                                     ref_gapped_seq, ref_genome_id, ref_seq_record, local_genome_names):
    """
    Compute similarity vectors and conservation percentages for the given gapped positions,
    based on the gapped motif sequence.
    """
    
    
    vectors = {}
    aligned_sequences = []
    total_conservation = 0
    genomes_with_data = 0
    motif_length = len(ref_gapped_seq.replace('-', ''))  # True length of the motif (excluding gaps)

    # Compute ungapped positions before motif start in reference genome
    ref_seq_gapped = str(ref_seq_record.seq)
    ref_strand = ref_seq_record.annotations.get('strand', '+')
    ref_start = int(ref_seq_record.annotations.get('start', 0))
    ref_size = int(ref_seq_record.annotations.get('size', 0))
    
    if ref_strand == -1:
        ref_strand = '-'
    elif ref_strand == 1:
        ref_strand = '+'

    # Compute ungapped positions before motif start
    ungapped_positions_before_motif = len([c for c in ref_seq_gapped[:gapped_start] if c != '-'])
    motif_ungapped_length = len(ref_gapped_seq.replace('-', ''))

    genomic_start = ref_start + ungapped_positions_before_motif + 1
    genomic_end = genomic_start + motif_ungapped_length - 1
    # Compute genomic start and end positions for reference genome
    # if ref_strand == '+':
    #     genomic_start = ref_start + ungapped_positions_before_motif + 1
    #     genomic_end = genomic_start + motif_ungapped_length - 1
    # else:
    #     genomic_end = ref_start + ref_size - ungapped_positions_before_motif
    #     genomic_start = genomic_end - motif_ungapped_length + 1

    # Extract chromosome information
    ref_chromosome = '.'.join(ref_seq_record.id.split('.')[1:])  # Skip the genome ID

    for seq_record in block:
        genome_id = seq_record.id.split('.')[0]
        if genome_id == ref_genome_id:
            continue
        local_genome_names.add(genome_id)  # Collect genome IDs
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
            seq_id_parts = seq_record.id.split('.')
            seq_chromosome = '.'.join(seq_id_parts[1:])  # Skip the genome ID

            if seq_strand == -1:
                seq_strand = '-'
            elif seq_strand == 1:
                seq_strand = '+'

            # Compute ungapped positions before motif start in this genome
            seq_gapped_full = str(seq_record.seq)
            ungapped_positions_before_motif_seq = len([c for c in seq_gapped_full[:gapped_start] if c != '-'])
            motif_ungapped_length_seq = len(seq_gapped_fragment.replace('-', ''))

            seq_genomic_start = seq_start + ungapped_positions_before_motif_seq + 1
            seq_genomic_end = seq_genomic_start + motif_ungapped_length_seq - 1

            # if seq_strand == '+':
            #     seq_genomic_start = seq_start + ungapped_positions_before_motif_seq + 1
            #     seq_genomic_end = seq_genomic_start + motif_ungapped_length_seq - 1
            # else:
            #     seq_genomic_end = seq_start + seq_size - ungapped_positions_before_motif_seq
            #     seq_genomic_start = seq_genomic_end - motif_ungapped_length_seq + 1

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

def convert_numpy_types(obj):
    """Recursively convert NumPy data types to native Python types."""
    if isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(v) for v in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

def search_patterns_in_block(block, search_type, A, patterns, genome_ids, block_no, search_in, genome_files, reverse_complement, process_id, local_genome_names, output_identifier):
    """Search for motifs in the sequences based on the specified search type."""
    try:
        # Determine the reference genome in this block
        ref_seq_record = block[0]  # Reference genome is always the first sequence
        local_genome_names.add(ref_seq_record.id.split('.')[0])

        if search_in == 'reference':
            # Search in the reference genome only
            sequences_to_search = [ref_seq_record]
            logging.debug("Search type: Reference genome only.")
        else:
            # Search in the specified genome
            genome_to_search = search_in
            sequences_to_search = [seq for seq in block if seq.id.split('.')[0] == genome_to_search]
            if not sequences_to_search:
                logging.debug(f"No sequences found for genome {genome_to_search} in block {block_no}.")
            logging.debug(f"Search type: Genome {genome_to_search}.")

        ref_genome_id = ref_seq_record.id.split('.')[0]

        for seq_record in sequences_to_search:
            genome_id = seq_record.id.split('.')[0]
            local_genome_names.add(genome_id)  # Collect genome IDs

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
            seq_ungapped_orig = seq_gapped_orig.replace('-', '').upper()

            if search_type == 'pwm':
                # PWM search
                for pwm_data in patterns:
                    pssm = pwm_data['pwm']
                    pwm_name = pwm_data['name']
                    threshold = pwm_data['threshold']
                    motif_identifier = pwm_data['identifier']
                    is_reverse_complement = pwm_data['is_reverse_complement']

                    m = pssm.length
                    # Ensure the sequence length is sufficient
                    if len(seq_ungapped_orig) < m:
                        continue

                    # Search for hits using precomputed threshold
                    hits = pssm.search(seq_ungapped_orig, threshold=threshold , both=False)
                    hits = list(hits)

                    for position, score in hits:
                        pos = position
                        # Get the matched sequence
                        ungapped_sequence = seq_ungapped_orig[pos:pos + m]
                        # Map ungapped positions to gapped positions
                        ungapped_to_gapped = [i for i, c in enumerate(seq_gapped_orig) if c != '-']
                        if pos + m - 1 >= len(ungapped_to_gapped):
                            continue
                        gapped_start = ungapped_to_gapped[pos]
                        gapped_end = ungapped_to_gapped[pos + m - 1] + 1

                        # Extract gapped sequence
                        gapped_sequence = seq_gapped_orig[gapped_start:gapped_end]
                        # Determine motif strand
                        if is_reverse_complement:
                            motif_strand = '-' if sequence_strand == '+' else '+'
                            hit_id_suffix = '_rc'
                        else:
                            motif_strand = sequence_strand
                            hit_id_suffix = ''

                        # Extract reference gapped sequence
                        ref_seq_gapped = str(ref_seq_record.seq)
                        ref_gapped_seq = ref_seq_gapped[gapped_start:gapped_end]

                        # Compute conservation and genomic coordinates
                        vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                            block, genome_ids, gapped_start, gapped_end, ref_gapped_seq, ref_genome_id, ref_seq_record, local_genome_names
                        )
                        # Collect motif hit data
                        motif_hit_data = {
                            "hit_id": f"pwm_{genome_id}_{motif_identifier}{hit_id_suffix}_{block_no}_{pos}",
                            "block_no": block_no,
                            "reference_genome": ref_genome_id,
                            "chromosome": ref_chromosome,
                            "genomic_start": genomic_start,
                            "genomic_end": genomic_end,
                            "strand": sequence_strand,
                            "motif_type": "pwm",
                            "motif_name": motif_identifier,
                            "motif_length": m,
                            "gapped_motif": gapped_sequence,
                            "gapped_start": gapped_start + 1,
                            "gapped_end": gapped_end,
                            "motif": ungapped_sequence,
                            "ungapped_start": pos + 1,
                            "ungapped_end": pos + m,
                            "motif_strand": motif_strand,
                            "score": float(score),
                            "conservation": {
                                "value": float(conservation_value),
                                "other_genomes": aligned_sequences
                            }
                        }

                        # Convert NumPy types to native Python types
                        motif_hit_data_converted = convert_numpy_types(motif_hit_data)

                        # Create output files for genomes and PWMs if not already created
                        key = (genome_id, motif_identifier)
                        if key not in genome_files:
                            genome_dir = os.path.join('results', genome_id)
                            os.makedirs(genome_dir, exist_ok=True)
                            output_file = os.path.join(genome_dir, f'{output_identifier}_{motif_identifier}_motif_hits_process_{process_id}.json')
                            try:
                                genome_file = open(output_file, 'a')
                                genome_files[key] = genome_file
                                logging.info(f"Process {process_id}: Created output file for genome {genome_id}, PWM {motif_identifier}")
                            except Exception as e:
                                logging.error(f"Process {process_id}: Error creating file {output_file}: {e}", exc_info=True)
                                continue

                        genome_file = genome_files[key]  # Get the file handle for the genome and motif

                        # Write the motif hit data to the JSON file
                        try:
                            json.dump(motif_hit_data_converted, genome_file)
                            genome_file.write('\n')  # Write newline to separate JSON objects
                        except Exception as e:
                            logging.error(f"Error writing motif hit: {e}", exc_info=True)

            elif search_type == 'kmer':
                # K-mer search using Aho-Corasick algorithm
                for end_index, (_, kmer_data) in A.iter(seq_ungapped_orig):
                    start_index = end_index - len(kmer_data['kmer']) + 1
                    pos = start_index
                    kmer = kmer_data['kmer']
                    original_kmer = kmer_data['original_kmer']
                    is_reverse_complement = kmer_data['is_reverse_complement']

                    kmer_len = len(kmer)
                    ungapped_sequence = seq_ungapped_orig[pos:pos + kmer_len]

                    # Map ungapped positions to gapped positions
                    ungapped_to_gapped = [i for i, c in enumerate(seq_gapped_orig) if c != '-']
                    if pos + kmer_len - 1 >= len(ungapped_to_gapped):
                        continue
                    gapped_start = ungapped_to_gapped[pos]
                    gapped_end = ungapped_to_gapped[pos + kmer_len - 1] + 1

                    # Extract gapped sequence
                    gapped_sequence = seq_gapped_orig[gapped_start:gapped_end]

                    # Determine motif strand
                    if is_reverse_complement:
                        motif_strand = '-' if sequence_strand == '+' else '+'
                        hit_id_suffix = '_rc'
                    else:
                        motif_strand = sequence_strand
                        hit_id_suffix = ''

                    # Extract reference gapped sequence
                    ref_seq_gapped = str(ref_seq_record.seq)
                    ref_gapped_seq = ref_seq_gapped[gapped_start:gapped_end]

                    # Compute conservation and genomic coordinates
                    vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                        block, genome_ids, gapped_start, gapped_end, ref_gapped_seq, ref_genome_id, ref_seq_record, local_genome_names
                    )

                    # Collect motif hit data
                    motif_hit_data = {
                        "hit_id": f"kmer_{genome_id}_{original_kmer}{hit_id_suffix}_{block_no}_{pos}",
                        "block_no": block_no,
                        "reference_genome": ref_genome_id,
                        "chromosome": ref_chromosome,
                        "genomic_start": genomic_start,
                        "genomic_end": genomic_end,
                        "strand": sequence_strand,
                        "motif_type": "kmer",
                        "motif_name": original_kmer,
                        "motif_length": kmer_len,
                        "gapped_motif": gapped_sequence,
                        "gapped_start": gapped_start + 1,
                        "gapped_end": gapped_end,
                        "motif": ungapped_sequence,
                        "ungapped_start": pos + 1,
                        "ungapped_end": pos + kmer_len,
                        "motif_strand": motif_strand,
                        "conservation": {
                            "value": float(conservation_value),
                            "other_genomes": aligned_sequences
                        }
                    }

                    # Convert NumPy types to native Python types
                    motif_hit_data_converted = convert_numpy_types(motif_hit_data)

                    # Create output files for genomes and K-mers if not already created
                    key = genome_id
                    if key not in genome_files:
                        genome_dir = os.path.join('results', genome_id)
                        os.makedirs(genome_dir, exist_ok=True)
                        output_file = os.path.join(genome_dir, f'{output_identifier}_motif_hits_process_{process_id}.json')
                        try:
                            genome_file = open(output_file, 'a')
                            genome_files[key] = genome_file
                            logging.info(f"Process {process_id}: Created output file for genome {genome_id}")
                        except Exception as e:
                            logging.error(f"Process {process_id}: Error creating file {output_file}: {e}", exc_info=True)
                            continue

                    genome_file = genome_files[key]  # Get the file handle for the genome

                    # Write the motif hit data to the JSON file
                    try:
                        json.dump(motif_hit_data_converted, genome_file)
                        genome_file.write('\n')  # Write newline to separate JSON objects
                    except Exception as e:
                        logging.error(f"Error writing motif hit: {e}", exc_info=True)

            elif search_type == 'regex':
                # Regex search
                for regex_data in patterns:
                    pattern = regex_data['pattern']
                    regex_name = regex_data['regex']
                    regex_identifier = regex_data['identifier']
                    is_reverse_complement = regex_data['is_reverse_complement']

                    if is_reverse_complement:
                        # seq_gapped_orig = str(Seq(seq_gapped_orig).reverse_complement())
                        sequence_to_search = str(Seq(seq_ungapped_orig).reverse_complement())
                        total_len = len(seq_ungapped_orig)
                        motif_strand = '-' if sequence_strand == '+' else '+'
                        hit_id_suffix = '_rc'
                    else:
                        sequence_to_search = seq_ungapped_orig
                        motif_strand = sequence_strand
                        hit_id_suffix = ''

                    # Search for matches
                    for match in pattern.finditer(sequence_to_search):
                        start_index = match.start()
                        end_index = match.end()
                        match_len = end_index - start_index
                        ungapped_sequence = sequence_to_search[start_index:end_index]

                        if is_reverse_complement:
                            pos_in_original = total_len - end_index
                        else:
                            pos_in_original = start_index

                        # Map ungapped positions to gapped positions
                        ungapped_to_gapped = [i for i, c in enumerate(seq_gapped_orig) if c != '-']
                        if pos_in_original + match_len - 1 >= len(ungapped_to_gapped):
                            continue
                        gapped_start = ungapped_to_gapped[pos_in_original]
                        gapped_end = ungapped_to_gapped[pos_in_original + match_len - 1] + 1

                        # Extract gapped sequence
                        gapped_sequence = seq_gapped_orig[gapped_start:gapped_end]

                        # Extract reference gapped sequence
                        ref_seq_gapped = str(ref_seq_record.seq)
                        ref_gapped_seq = ref_seq_gapped[gapped_start:gapped_end]

                        # Compute conservation and genomic coordinates
                        vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                            block, genome_ids, gapped_start, gapped_end, ref_gapped_seq, ref_genome_id, ref_seq_record, local_genome_names
                        )

                        # Collect motif hit data
                        motif_hit_data = {
                            "hit_id": f"regex_{genome_id}_{regex_identifier}{hit_id_suffix}_{block_no}_{pos_in_original}",
                            "block_no": block_no,
                            "reference_genome": ref_genome_id,
                            "chromosome": ref_chromosome,
                            "genomic_start": genomic_start,
                            "genomic_end": genomic_end,
                            "strand": sequence_strand,
                            "motif_type": "regex",
                            "motif_name": regex_identifier,
                            "motif_length": match_len,
                            "gapped_motif": gapped_sequence,
                            "gapped_start": gapped_start + 1,
                            "gapped_end": gapped_end,
                            "ungapped_start": pos_in_original + 1,
                            "ungapped_end": pos_in_original + match_len,
                            "motif": ungapped_sequence,
                            "motif_strand": motif_strand,
                            "conservation": {
                                "value": float(conservation_value),
                                "other_genomes": aligned_sequences
                            }
                        }

                        # Convert NumPy types to native Python types
                        motif_hit_data_converted = convert_numpy_types(motif_hit_data)

                        # Create output files for genomes and regex patterns if not already created
                        key = (genome_id, regex_identifier)
                        if key not in genome_files:
                            genome_dir = os.path.join('results', genome_id)
                            os.makedirs(genome_dir, exist_ok=True)
                            output_file = os.path.join(genome_dir, f'{output_identifier}_{regex_identifier}_motif_hits_process_{process_id}.json')
                            try:
                                genome_file = open(output_file, 'a')
                                genome_files[key] = genome_file
                                logging.info(f"Process {process_id}: Created output file for genome {genome_id}, regex {regex_identifier}")
                            except Exception as e:
                                logging.error(f"Process {process_id}: Error creating file {output_file}: {e}", exc_info=True)
                                continue

                        genome_file = genome_files[key]  # Get the file handle for the genome and regex

                        # Write the motif hit data to the JSON file
                        try:
                            json.dump(motif_hit_data_converted, genome_file)
                            genome_file.write('\n')  # Write newline to separate JSON objects
                        except Exception as e:
                            logging.error(f"Error writing motif hit: {e}", exc_info=True)

            else:
                logging.error(f"Unknown search type: {search_type}")
                sys.exit(1)

    except Exception as e:
        logging.error(f"Error processing block {block_no}: {e}", exc_info=True)
        sys.exit(1)  # Stop the script

def write_motif_hit(genome_file, motif_hit_data):
    """Write a motif hit to the genome's JSON file."""
    try:
        json.dump(motif_hit_data, genome_file)
        genome_file.write('\n')  # Write newline to separate JSON objects
    except Exception as e:
        logging.error(f"Error writing motif hit: {e}", exc_info=True)

def _initialize_automaton(patterns):
    A = ahocorasick.Automaton()  # Updated to use 'ahocorasick'
    for idx, kmer_data in enumerate(patterns):
        kmer = kmer_data['kmer']
        kmer_data['index'] = idx
        A.add_word(kmer, (idx, kmer_data))
    A.make_automaton()
    return A

def process_file_chunk(maf_file, start_pos, end_pos, search_type, patterns, genome_ids, search_in, reverse_complement, process_id, unique_genome_names, output_identifier):
    """Process a chunk of the MAF file."""
    genome_files = {}
    block_no = 0  # Initialize block number
    local_genome_names = set()
    if search_type == "kmer":
        A = _initialize_automaton(patterns)
    else:
        A = None

    try:
        with open(maf_file, 'rb') as handle:
            handle.seek(start_pos)
            current_position = handle.tell()
            line = handle.readline()
            while current_position < end_pos and line:
                if not line.startswith(b'a'):
                    current_position = handle.tell()
                    line = handle.readline()
                    continue
                # Start of a block
                block_data = line.decode('utf-8')
                current_position = handle.tell()
                line = handle.readline()
                while line and not line.startswith(b'a') and current_position < end_pos:
                    block_data += line.decode('utf-8')
                    current_position = handle.tell()
                    line = handle.readline()
                # Parse the block
                block = parse_block_from_string(block_data)
                if block:
                    # Search patterns in the block
                    search_patterns_in_block(block, search_type, A, patterns, genome_ids, block_no, search_in, genome_files, reverse_complement, process_id, local_genome_names, output_identifier)
                block_no += 1
    except Exception as e:
        logging.error(f"Error in process {process_id}: {e}", exc_info=True)
    finally:
        # Close all genome files
        for genome_file in genome_files.values():
            try:
                genome_file.close()
            except Exception as e:
                logging.error(f"Error closing genome file: {e}", exc_info=True)
        # Extend the unique genome names
        unique_genome_names.extend(local_genome_names)

def merge_results(results_dir, output_identifier):
    """Merge results from all processes into per-genome JSON files."""
    logging.info("Merging results from all processes.")
    # Get list of genome directories
    genome_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    for genome_id in genome_dirs:
        genome_dir = os.path.join(results_dir, genome_id)
        # Collect all process-specific files in the genome directory
        process_files = [f for f in os.listdir(genome_dir) if f.endswith('.json')]
        files_by_motif = {}
        for f in process_files:
            match = re.match(rf'{output_identifier}_(.+?)_motif_hits_process_\d+\.json', f)
            if match:
                motif_identifier = match.group(1)
                files_by_motif.setdefault(motif_identifier, []).append(f)
            else:
                # Handle k-mer files
                if f.startswith(f'{output_identifier}_motif_hits_process_'):
                    files_by_motif.setdefault(output_identifier, []).append(f)

        for motif_identifier, files in files_by_motif.items():
            merged_data = []
            files.sort(key=lambda x: int(re.findall(r'\d+', x)[-1]) if re.findall(r'\d+', x) else 0)
            for motif_file in files:
                motif_file_path = os.path.join(genome_dir, motif_file)
                if os.path.exists(motif_file_path):
                    with open(motif_file_path, 'r') as f:
                        for line in f:
                            try:
                                data = json.loads(line)
                                merged_data.append(data)
                            except json.JSONDecodeError:
                                logging.error(f"Error decoding JSON from {motif_file_path}")
            if merged_data:
                # Write merged data to final result file
                final_result_file = os.path.join(genome_dir, f'{output_identifier}_{motif_identifier}_motif_hits.json')
                try:
                    with open(final_result_file, 'w') as f:
                        json.dump(merged_data, f, indent=4)
                    logging.info(f"Merged results for genome {genome_id}, motif {motif_identifier} into {final_result_file}")
                except Exception as e:
                    logging.error(f"Error writing merged results for genome {genome_id}, motif {motif_identifier}: {e}", exc_info=True)
                # Optionally, remove the process-specific files
                for motif_file in files:
                    try:
                        os.remove(os.path.join(genome_dir, motif_file))
                        logging.info(f"Removed temporary file {motif_file} for genome {genome_id}, motif {motif_identifier}")
                    except Exception as e:
                        logging.error(f"Error removing temporary file {motif_file}: {e}", exc_info=True)
            else:
                logging.info(f"No data found for genome {genome_id}, motif {motif_identifier}")

def generate_csv(results_dir, unique_genome_names, output_identifier):
    """Generate CSV files with motif hits and conservation data."""
    logging.info("Generating CSV files.")
    unique_genome_names = set(unique_genome_names)

    # Get list of genome directories
    genome_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]

    for genome_id in genome_dirs:
        genome_dir = os.path.join(results_dir, genome_id)
        # Collect all merged JSON files in the genome directory
        merged_files = [f for f in os.listdir(genome_dir) if f.endswith('_motif_hits.json')]
        for merged_file in merged_files:
            motif_identifier_match = re.match(rf'{output_identifier}_(.+?)_motif_hits\.json', merged_file)
            if motif_identifier_match:
                motif_identifier = motif_identifier_match.group(1)
            else:
                motif_identifier = output_identifier  # For k-mer files

            motif_file = os.path.join(genome_dir, merged_file)
            all_motif_hits = {}
            if os.path.exists(motif_file):
                with open(motif_file, 'r') as f:
                    try:
                        data = json.load(f)
                        for motif_hit in data:
                            hit_id = motif_hit['hit_id']
                            all_motif_hits[hit_id] = motif_hit
                    except json.JSONDecodeError:
                        logging.error(f"Error decoding JSON from {motif_file}")
                    except Exception as e:
                        logging.error(f"Error reading {motif_file}: {e}", exc_info=True)
            else:
                logging.info(f"No motif_hits.json found for genome {genome_id}, motif {motif_identifier}")

            # Now, we have all motif hits in all_motif_hits
            # Build the CSV
            csv_file = os.path.join(genome_dir, f'{output_identifier}_{motif_identifier}_motif_hits.csv')
            fieldnames = ['motif_hit_info', 'motif_name', 'motif_type', 'motif_strand', 'score'] + list(unique_genome_names)
            try:
                with open(csv_file, 'w', newline='') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()

                    for hit_id, motif_hit in all_motif_hits.items():
                        # First columns: motif_hit_info and details
                        motif_hit_info = f"{motif_hit['chromosome']}:{motif_hit['genomic_start']}-{motif_hit['genomic_end']},{motif_hit['motif']}"
                        row = {
                            'motif_hit_info': motif_hit_info,
                            'motif_name': motif_hit['motif_name'],
                            'motif_type': motif_hit['motif_type'],
                            'motif_strand': motif_hit['motif_strand'],
                            'score': motif_hit.get('score')
                        }
                        # Initialize all genome columns to None
                        for genome_name in unique_genome_names:
                            row[genome_name] = None
                        # Fill in the data
                        conservation_data = motif_hit.get('conservation', {})
                        other_genomes = conservation_data.get('other_genomes', [])
                        for genome_data in other_genomes:
                            genome_name = genome_data['genome_id']
                            if genome_name in unique_genome_names:
                                vector = genome_data.get('vector')
                                conservation_pct = genome_data.get('conservation')
                                row[genome_name] = f"{vector},{conservation_pct}"
                        # Write the row
                        writer.writerow(row)
                logging.info(f"CSV file generated at {csv_file}")
            except Exception as e:
                logging.error(f"Error writing CSV file: {e}", exc_info=True)

    logging.info("CSV generation completed.")

def purge_directory(directory):
    if os.path.exists(directory):
        logging.info(f"Purging {directory}...")
        shutil.rmtree(directory)
        logging.info(f"{directory} has been purged.")
    else:
        logging.info(f"{directory} does not exist.")

def main():

    parser = argparse.ArgumentParser(description="Search uncompressed MAF files for specific motifs using regex patterns, K-mers, or PWMs.")

    parser.add_argument('--maf_file', required=True, help="Path to the uncompressed MAF file.")
    parser.add_argument('--genome_ids', help="Path to the genome IDs file (optional).")
    parser.add_argument('--search_in', default='reference',
                        help="Specify the genome to search motifs in (default: 'reference').")
    parser.add_argument('--reverse_complement', choices=['yes', 'no'], default='no',
                        help="Specify whether to search motifs on both strands.")
    parser.add_argument('--pvalue_threshold', type=float, default=1e-4,
                        help="P-value threshold for PWM matches. Default is 1e-4.")
    parser.add_argument('--processes', type=int, default=1,
                        help="Number of parallel processes to use. Default is 1.")
    parser.add_argument('--background_frequencies', nargs=4, type=float,
                    metavar=('A_FREQ', 'C_FREQ', 'G_FREQ', 'T_FREQ'),
                    help="Background nucleotide frequencies for A, C, G, T. They must sum to 1.")
    parser.add_argument('--purge_results_dir', action='store_true', help='Purge the results directory before running the script.')

    # Mutually exclusive group for search types
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--regexes', nargs='+', help="Regex patterns to search for in the MAF file.")
    group.add_argument('--kmers', help="Path to the K-mers file (one per line).")
    group.add_argument('--jaspar_file', help="Path to the PWM file in JASPAR format.")

    args = parser.parse_args()


    logging.info("Starting processing.")

    maf_file = args.maf_file
    genome_ids_file = args.genome_ids
    search_in = args.search_in
    reverse_complement = args.reverse_complement == 'yes'
    pvalue_threshold = args.pvalue_threshold
    num_processes = args.processes

    results_dir = os.path.join(os.getcwd(), 'results')
    
    if args.purge_results_dir:
        logging.info("Purging the results directory...")
        purge_directory(results_dir)
    
    # Add your code here to purge the directory

    os.makedirs(results_dir, exist_ok=True)

    genome_ids = load_genome_ids(genome_ids_file)

    # Process background frequencies
    background_frequencies = None
    if args.background_frequencies:
        total = sum(args.background_frequencies)
        if not np.isclose(total, 1.0):
            logging.error("The background frequencies must sum to 1.")
            sys.exit(1)
        nucleotides = ['A', 'C', 'G', 'T']
        background_frequencies = dict(zip(nucleotides, args.background_frequencies))
        logging.info(f"Using background frequencies from command-line argument: {background_frequencies}")
    else:
        background_frequencies = None

    # Determine search type and load patterns
    search_type = None
    patterns = None
    output_identifier = None

    if args.regexes:
        search_type = 'regex'
        regex_patterns = []
        output_identifier = 'regex_search'
        for regex in args.regexes:
            try:
                pattern = re.compile(regex, re.IGNORECASE)
                regex_identifier = re.sub(r'\W+', '_', regex)
                regex_patterns.append({'pattern': pattern, 'regex': regex,
                                       'identifier': regex_identifier,
                                       'is_reverse_complement': False})
                if reverse_complement:
                    regex_patterns.append({'pattern': pattern, 'regex': regex,
                                           'identifier': regex_identifier,
                                           'is_reverse_complement': True})
            except re.error as e:
                logging.error(f"Invalid regex pattern '{regex}': {e}")
                sys.exit(1)
        patterns = regex_patterns
        logging.info(f"Loaded {len(patterns)} regex patterns.")
    elif args.kmers:
        search_type = 'kmer'
        kmer_patterns = []
        kmers_file_name = os.path.splitext(os.path.basename(args.kmers))[0]
        output_identifier = kmers_file_name
        try:
            with open(args.kmers, 'r') as f:
                kmers = list(set([line.strip().upper() for line in f if line.strip()]))
            for kmer in kmers:
                kmer_patterns.append({'kmer': kmer, 'identifier': kmer,
                                      'original_kmer': kmer, 'is_reverse_complement': False})
                if reverse_complement:
                    rc_kmer = str(Seq(kmer).reverse_complement())
                    kmer_patterns.append({'kmer': rc_kmer, 'identifier': kmer,
                                          'original_kmer': kmer, 'is_reverse_complement': True})
            kmer_patterns = list({(d['kmer'], d['is_reverse_complement']): d for d in kmer_patterns}.values())
            patterns = kmer_patterns
            logging.info(f"Loaded {len(patterns)} K-mers (including reverse complements).")
        except Exception as e:
            logging.error(f"Error loading K-mers from {args.kmers}: {e}", exc_info=True)
            sys.exit(1)
    elif args.jaspar_file:
        search_type = 'pwm'
        pwm_motifs = []
        jaspar_file_name = os.path.splitext(os.path.basename(args.jaspar_file))[0]
        output_identifier = jaspar_file_name
        try:
            with open(args.jaspar_file, 'r') as f:
                for motif in motifs.parse(f, 'jaspar'):
                    motif_identifier = motif.matrix_id.strip()
                    pwm = motif.counts.normalize(pseudocounts=0.1)
                    pssm = pwm.log_odds()

                    # Set background frequencies
                    if background_frequencies:
                        pssm.background = background_frequencies
                        logging.info(f"Using background frequencies from command-line argument for motif {motif_identifier}: {background_frequencies}")
                    elif motif.background:
                        pssm.background = motif.background
                        logging.info(f"Using background frequencies from JASPAR file for motif {motif_identifier}: {motif.background}")
                    else:
                        pssm.background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
                        logging.info(f"Using uniform background frequencies for motif {motif_identifier}")

                    threshold_score = compute_threshold_score(pssm, pvalue_threshold)
                    pwm_motifs.append({'pwm': pssm, 'name': motif.name,
                                       'threshold': threshold_score,
                                       'identifier': motif_identifier,
                                       'is_reverse_complement': False})
                    if reverse_complement:
                        pssm_rc = pssm.reverse_complement()
                        pssm_rc.background = pssm.background  # Ensure background frequencies are set
                        threshold_rc = compute_threshold_score(pssm_rc, pvalue_threshold)
                        pwm_motifs.append({'pwm': pssm_rc, 'name': motif.name,
                                           'threshold': threshold_rc,
                                           'identifier': motif_identifier,
                                           'is_reverse_complement': True})
            patterns = pwm_motifs
            
            logging.info(f"Loaded {len(patterns)} PWMs (including reverse complements).")
        except Exception as e:
            logging.error(f"Error loading PWMs from {args.jaspar_file}: {e}", exc_info=True)
            sys.exit(1)
    else:
        logging.error("Error: One of --regexes, --kmers, or --jaspar_file must be provided.")
        sys.exit(1)

    chunk_ranges = assign_file_chunks(maf_file, num_processes)

    manager = multiprocessing.Manager()
    unique_genome_names = manager.list()
    processes = []
    for i, (start_pos, end_pos) in enumerate(chunk_ranges):
        p = multiprocessing.Process(target=process_file_chunk, args=(
            maf_file, start_pos, end_pos, search_type, patterns,
            genome_ids, search_in, reverse_complement, i,
            unique_genome_names, output_identifier))
        processes.append(p)
        p.start()
        logging.info(f"Started process {i} with PID {p.pid}")

    for i, p in enumerate(processes):
        p.join()
        logging.info(f"Process {i} with PID {p.pid} has completed.")

    unique_genome_names = list(set(unique_genome_names))
    logging.info(f"Unique genome names encountered: {unique_genome_names}")

    # Merge results from all processes
    merge_results(results_dir, output_identifier)

    # Generate CSV files
    generate_csv(results_dir, unique_genome_names, output_identifier)

    logging.info("Completed processing.")

if __name__ == "__main__":
    main()

