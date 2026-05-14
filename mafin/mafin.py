#!/usr/bin/env python3

import sys
import re
import argparse
import json
import logging
import os
import multiprocessing
from Bio import motifs
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq
import random
import numpy as np
import csv
import shutil
import ahocorasick

# Try to use a faster JSON library if available
try:
    import ujson as fast_json
except ImportError:
    try:
        import orjson as fast_json
    except ImportError:
        fast_json = json

def fast_dumps(data):
    """Fallback-safe high speed JSON dumping"""
    encoded = fast_json.dumps(data)
    return encoded.decode('utf-8') if isinstance(encoded, bytes) else encoded

# We'll configure logging dynamically based on --verbose
logger = logging.getLogger(__name__)

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)


def load_genome_ids(genome_ids_file):
    """Load genome IDs from a given file."""
    if genome_ids_file:
        logger.info(f"Loading genome IDs from {genome_ids_file}")
        try:
            with open(genome_ids_file, 'r') as file:
                # BUD Optimization: Use a Hash Set instead of an Array for O(1) lookups
                genome_ids = {line.strip() for line in file if line.strip()}
            logger.info(f"Loaded {len(genome_ids)} genome IDs.")
            return genome_ids
        except Exception as e:
            logger.error(f"Error loading genome IDs from {genome_ids_file}: {e}", exc_info=True)
            sys.exit(1)
    return None  # If not provided, we consider all genomes


def assign_file_chunks(maf_file, num_processes):
    """
    Assign chunks of the MAF file to processes based on block boundaries.
    Each chunk ends just before the next line that starts with 'a'.
    """
    file_size = os.path.getsize(maf_file)
    # Edge case: if num_processes is 0 or file_size is 0
    if num_processes < 1 or file_size == 0:
        return [(0, file_size)]

    chunk_size = file_size // num_processes
    chunk_positions = [0]  # The first boundary is always 0

    with open(maf_file, 'rb') as f:
        offset = 0
        for _ in range(1, num_processes):
            nominal_start = offset + chunk_size
            if nominal_start >= file_size:
                chunk_positions.append(file_size)
                offset = file_size
                continue

            f.seek(nominal_start)
            found_a = False
            while True:
                current_pos = f.tell()
                line = f.readline()
                if not line:
                    chunk_positions.append(file_size)
                    offset = file_size
                    found_a = True
                    break

                if line.startswith(b'a'):
                    chunk_positions.append(current_pos)
                    offset = current_pos
                    found_a = True
                    break

            if not found_a:
                break

        while len(chunk_positions) < num_processes + 1:
            if chunk_positions[-1] != file_size:
                chunk_positions.append(file_size)
            else:
                chunk_positions.append(file_size)

    chunk_ranges = []
    for i in range(num_processes):
        start_pos = chunk_positions[i]
        end_pos = chunk_positions[i + 1]
        logger.info(f"Assigned chunk {i}: Start at {start_pos}, End at {end_pos}")
        chunk_ranges.append((start_pos, end_pos))

    return chunk_ranges


class SimpleSeqRecord:
    __slots__ = ['id', 'genome_id', 'chromosome', 'annotations', 'seq', 'seq_np']

    def __init__(self):
        self.id = ""
        self.genome_id = ""
        self.chromosome = ""
        self.annotations = {'start': 0, 'size': 0, 'strand': 1, 'srcSize': 0}
        self.seq = ""
        self.seq_np = None

    def update(self, obj_id, start, size, strand, srcSize, seq):
        self.id = obj_id
        parts = obj_id.split('.', 1)
        self.genome_id = parts[0]
        self.chromosome = parts[1] if len(parts) > 1 else 'unknown'
        self.annotations['start'] = start
        self.annotations['size'] = size
        self.annotations['strand'] = strand
        self.annotations['srcSize'] = srcSize
        self.seq = seq
        self.seq_np = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)


def parse_block_from_string(block_data, record_pool):
    """Parse a MAF block from a string into a list of SimpleSeqRecord objects using an object pool."""
    num_records = 0
    try:
        lines = block_data.strip().split('\n')
        for line in lines:
            if line.startswith('s '):
                parts = line.split()
                if len(parts) >= 7:
                    src = parts[1]
                    start = int(parts[2])
                    size = int(parts[3])
                    strand_str = parts[4]
                    strand = 1 if strand_str == '+' else -1
                    srcSize = int(parts[5])
                    text = parts[6].upper()
                    
                    if num_records < len(record_pool):
                        record_pool[num_records].update(src, start, size, strand, srcSize, text)
                    else:
                        rec = SimpleSeqRecord()
                        rec.update(src, start, size, strand, srcSize, text)
                        record_pool.append(rec)
                    
                    num_records += 1
                    
        return record_pool[:num_records] if num_records else None
    except Exception as e:
        logger.error(f"Error parsing block: {e}", exc_info=True)
        return None


def compute_threshold_score(pssm, pvalue, num_samples=10000):
    """Compute the threshold score for a given p-value by sampling random sequences."""
    logger.info(f"Computing threshold score for p-value: {pvalue} with background frequencies: {pssm.background}")
    scores = []
    for _ in range(num_samples):
        seq = generate_random_sequence(pssm.length, pssm.background)
        score = pssm.calculate(seq)
        scores.append(score)
    scores.sort()
    index = int((1 - pvalue) * num_samples)
    threshold_score = scores[index] if index < len(scores) else scores[-1]
    logger.info(f"Threshold score corresponding to p-value {pvalue}: {threshold_score}")
    return threshold_score


def generate_random_sequence(length, background):
    """Generate a random sequence based on the background nucleotide frequencies."""
    nucleotides = ['A', 'C', 'G', 'T']
    frequencies = [background.get(nuc, 0.25) for nuc in nucleotides]
    total_freq = sum(frequencies)
    if not np.isclose(total_freq, 1.0):
        frequencies = [freq / total_freq for freq in frequencies]
    return Seq(''.join(random.choices(nucleotides, frequencies, k=length)))


def compute_vectors_and_conservation(block,
                                     ref_seq_record,
                                     genome_ids,
                                     gapped_start,
                                     gapped_end,
                                     ref_genome_id,
                                     local_genome_names,
                                     ref_ungapped_start,
                                     ref_ungapped_length):
    """
    Compute similarity vectors and conservation percentages for the given gapped positions.
    """

    vectors = {}
    aligned_sequences = []
    total_conservation = 0.0
    genomes_with_data = 0

    # Locate reference sequence metadata
    ref_start = int(ref_seq_record.annotations.get('start', 0))

    genomic_start = int(ref_start + ref_ungapped_start)
    genomic_end = int(genomic_start + ref_ungapped_length - 1)
    ref_chromosome = ref_seq_record.chromosome

    # Extract reference sequence byte-array directly from pre-computed view
    r_arr_full = ref_seq_record.seq_np[gapped_start:gapped_end]
    r_len = len(r_arr_full)

    for meta in block:
        genome_id = meta.genome_id
        if genome_id == ref_genome_id:
            continue
        local_genome_names.add(genome_id)

        # If a list of genome_ids was specified, only process those
        if genome_ids is not None and genome_id not in genome_ids:
            continue

        # Zero-copy slicing of the pre-encoded numpy array
        o_arr = meta.seq_np[gapped_start:gapped_end]

        if len(o_arr) != r_len:
            # Pad with '-' (45) if it's shorter
            padded = np.full(r_len, 45, dtype=np.uint8)
            padded[:len(o_arr)] = o_arr
            o_arr = padded

        # Vectorized comparison logic utilizing AVX/SIMD CPU instructions and cache-aligned buffers
        valid_mask = ~((r_arr_full == 45) & (o_arr == 45))
        r_valid = r_arr_full[valid_mask]
        o_valid = o_arr[valid_mask]

        vector_length = len(r_valid)
        if vector_length == 0:
            continue

        vector_arr = np.full(vector_length, 48, dtype=np.uint8) # Fill with '0'
        
        matches_mask = (r_valid == o_valid)
        vector_arr[matches_mask] = 49 # '1'
        vector_arr[r_valid == 45] = 45 # '-'

        matches = np.count_nonzero(matches_mask)
        vector_str = vector_arr.tobytes().decode('ascii')
        conservation_pct = (matches / vector_length) * 100.0

        total_conservation += conservation_pct
        genomes_with_data += 1

        seq_gapped = meta.seq
        seq_gapped_fragment = seq_gapped[gapped_start:gapped_end]

        # O(1) Data mapping using precomputed cached coordinate arrays
        ungapped_positions_before_motif_seq = gapped_start - seq_gapped.count('-', 0, gapped_start)
        motif_ungapped_length_seq = len(seq_gapped_fragment) - seq_gapped_fragment.count('-')

        meta_start = int(meta.annotations.get('start', 0))
        seq_genomic_start = int(meta_start + ungapped_positions_before_motif_seq)
        seq_genomic_end = int(seq_genomic_start + motif_ungapped_length_seq - 1)
        
        mstrand = meta.annotations.get('strand', 1)
        mstrand_str = '-' if mstrand == -1 else '+'

        aligned_sequences.append({
            "genome_id": genome_id,
            "chromosome": meta.chromosome,
            "genomic_start": seq_genomic_start,
            "genomic_end": seq_genomic_end,
            "strand": mstrand_str,
            "aligned_sequence_gapped": seq_gapped_fragment,
            "vector": vector_str,
            "conservation": f"{conservation_pct:.2f}%",
            "gapped_start": gapped_start,
            "gapped_end": gapped_end - 1,
        })
        vectors[genome_id] = vector_str

    if genomes_with_data > 0:
        # Average conservation across all compared genomes
        avg_conservation_value = total_conservation / genomes_with_data
    else:
        avg_conservation_value = 0.0

    avg_conservation_value = round(avg_conservation_value, 2)

    return vectors, avg_conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome


def write_hits_to_disk(genome_dir, output_file, lines):
    """Background I/O task function to write motif hits."""
    try:
        os.makedirs(genome_dir, exist_ok=True)
        with open(output_file, 'a') as f:
            f.write('\n'.join(lines) + '\n')
    except Exception as e:
        logger.error(f"Error buffering motif hit to disk: {e}", exc_info=True)

def search_patterns_in_block(block,
                             search_type,
                             A,
                             patterns,
                             genome_ids,
                             block_no,
                             search_in,
                             genome_files,
                             reverse_complement,
                             process_id,
                             local_genome_names,
                             output_identifier,
                             io_executor):
    """Search for motifs in the sequences based on the specified search type."""
    try:
        
        if search_in == 'reference':
            ref_seq_record = block[0]
            local_genome_names.add(ref_seq_record.genome_id)
            sequences_to_search = [ref_seq_record]
        else:
            # If user gave a specific genome X, we find that in the block
            genome_to_search = search_in
            found_records = [seq for seq in block if seq.genome_id == genome_to_search]
            if not found_records:
                logger.debug(f"No sequences found for genome {genome_to_search} in block {block_no}.")
                return  # skip this block entirely

            ref_seq_record = found_records[0]
            ref_genome_id = ref_seq_record.genome_id
            local_genome_names.add(ref_genome_id)
            sequences_to_search = [ref_seq_record]

        ref_genome_id = ref_seq_record.genome_id

        for seq_record in sequences_to_search:
            genome_id = seq_record.genome_id
            local_genome_names.add(genome_id)

            seq_gapped = seq_record.seq
            seq_ungapped_orig = seq_gapped.replace('-', '')

            sstrand = seq_record.annotations.get('strand', 1)
            sequence_strand = '-' if sstrand == -1 else '+'

            ungapped_to_gapped = None 

            if search_type == 'pwm':
                # PWM search
                for pwm_data in patterns:
                    pssm = pwm_data['pwm']
                    threshold = pwm_data['threshold']
                    motif_identifier = pwm_data['identifier']
                    is_rc = pwm_data['is_reverse_complement']

                    m = pssm.length
                    if len(seq_ungapped_orig) < m:
                        continue

                    hits = pssm.search(seq_ungapped_orig, threshold=threshold, both=False)
                    hits = list(hits)

                    for position, score in hits:
                        score = round(float(score), 2)
                        pos = position
                        ungapped_sequence = seq_ungapped_orig[pos:pos + m]
                        
                        if ungapped_to_gapped is None:
                            ungapped_to_gapped = np.flatnonzero(seq_record.seq_np != 45)

                        if pos + m - 1 >= len(ungapped_to_gapped):
                            continue
                        # 0-based inclusive
                        gapped_start = int(ungapped_to_gapped[pos])
                        gapped_end = int(ungapped_to_gapped[pos + m - 1])

                        gapped_sequence = seq_gapped[gapped_start:gapped_end + 1]

                        if is_rc:
                            motif_strand = '-' if sequence_strand == '+' else '+'
                            hit_id_suffix = '_rc'
                        else:
                            motif_strand = sequence_strand
                            hit_id_suffix = ''

                        vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                            block,
                            ref_seq_record,
                            genome_ids,
                            gapped_start,
                            gapped_end + 1,
                            ref_genome_id,
                            local_genome_names,
                            pos,
                            m
                        )

                        motif_hit_data = {
                            "hit_id": f"pwm_{genome_id}_{motif_identifier}{hit_id_suffix}_{pos}",
                            "reference_genome": ref_genome_id,
                            "chromosome": ref_chromosome,
                            "genomic_start": genomic_start,
                            "genomic_end": genomic_end,
                            "strand": sequence_strand,
                            "motif_type": "pwm",
                            "motif_name": motif_identifier,
                            "motif_length": m,
                            "gapped_motif": gapped_sequence,
                            "gapped_start": gapped_start,
                            "gapped_end": gapped_end,
                            "motif": ungapped_sequence,
                            "ungapped_start": pos,
                            "ungapped_end": pos + m - 1,
                            "motif_strand": motif_strand,
                            "score": float(score),
                            "conservation": {
                                "value": float(conservation_value),
                                "other_genomes": aligned_sequences
                            }
                        }

                        key = (genome_id, motif_identifier)
                        if key not in genome_files:
                            genome_files[key] = []
                        
                        try:
                            genome_files[key].append(fast_dumps(motif_hit_data))
                            if len(genome_files[key]) >= 10000:
                                genome_dir = os.path.join('tmp', genome_id)
                                output_file = os.path.join(
                                    genome_dir,
                                    f'{output_identifier}_{motif_identifier}_motif_hits_process_{process_id}.tmp'
                                )
                                io_executor.submit(write_hits_to_disk, genome_dir, output_file, list(genome_files[key]))
                                genome_files[key].clear()
                        except Exception as e:
                            logger.error(f"Error buffering motif hit: {e}", exc_info=True)

            elif search_type == 'kmer':
                # K-mer search
                for end_index, (_, kmer_data) in A.iter(seq_ungapped_orig):
                    start_index = end_index - len(kmer_data['kmer']) + 1
                    pos = start_index
                    kmer = kmer_data['kmer']
                    original_kmer = kmer_data['original_kmer']
                    is_rc = kmer_data['is_reverse_complement']

                    kmer_len = len(kmer)
                    ungapped_sequence = seq_ungapped_orig[pos:pos + kmer_len]

                    if ungapped_to_gapped is None:
                        ungapped_to_gapped = np.flatnonzero(seq_record.seq_np != 45)

                    if pos + kmer_len - 1 >= len(ungapped_to_gapped):
                        continue

                    gapped_start = int(ungapped_to_gapped[pos])
                    gapped_end = int(ungapped_to_gapped[pos + kmer_len - 1])

                    gapped_sequence = seq_gapped[gapped_start:gapped_end + 1]

                    if is_rc:
                        motif_strand = '-' if sequence_strand == '+' else '+'
                        hit_id_suffix = '_rc'
                    else:
                        motif_strand = sequence_strand
                        hit_id_suffix = ''

                    vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                        block,
                        ref_seq_record,
                        genome_ids,
                        gapped_start,
                        gapped_end + 1,
                        ref_genome_id,
                        local_genome_names,
                        pos,
                        kmer_len
                    )

                    motif_hit_data = {
                        "hit_id": f"kmer_{genome_id}_{original_kmer}{hit_id_suffix}_{pos}",
                        "reference_genome": ref_genome_id,
                        "chromosome": ref_chromosome,
                        "genomic_start": genomic_start,
                        "genomic_end": genomic_end,
                        "strand": sequence_strand,
                        "motif_type": "kmer",
                        "motif_name": original_kmer,
                        "motif_length": kmer_len,
                        "gapped_motif": gapped_sequence,
                        "gapped_start": gapped_start,
                        "gapped_end": gapped_end,
                        "motif": ungapped_sequence,
                        "ungapped_start": pos,
                        "ungapped_end": pos + kmer_len - 1,
                        "motif_strand": motif_strand,
                        "conservation": {
                            "value": float(round(conservation_value, 2)),
                            "other_genomes": aligned_sequences
                        }
                    }

                    key = genome_id
                    if key not in genome_files:
                        genome_files[key] = []
                    
                    try:
                        genome_files[key].append(fast_dumps(motif_hit_data))
                        
                        if len(genome_files[key]) >= 10000:
                            genome_dir = os.path.join('tmp', genome_id)
                            output_file = os.path.join(
                                genome_dir,
                                f'{output_identifier}_motif_hits_process_{process_id}.tmp'
                            )
                            io_executor.submit(write_hits_to_disk, genome_dir, output_file, list(genome_files[key]))
                            genome_files[key].clear()
                    except Exception as e:
                        logger.error(f"Error buffering motif hit: {e}", exc_info=True)

            elif search_type == 'regex':
                # Regex search
                seq_ungapped_orig_rc = None
                for regex_data in patterns:
                    pattern = regex_data['pattern']
                    regex_identifier = regex_data['identifier']
                    is_rc = regex_data['is_reverse_complement']

                    if is_rc:
                        if seq_ungapped_orig_rc is None:
                            seq_ungapped_orig_rc = str(Seq(seq_ungapped_orig).reverse_complement())
                        sequence_to_search = seq_ungapped_orig_rc
                        total_len = len(seq_ungapped_orig)
                        motif_strand = '-' if sequence_strand == '+' else '+'
                        hit_id_suffix = '_rc'
                    else:
                        sequence_to_search = seq_ungapped_orig
                        motif_strand = sequence_strand
                        hit_id_suffix = ''

                    for match in pattern.finditer(sequence_to_search):
                        start_index = match.start()
                        end_index = match.end()
                        match_len = end_index - start_index
                        ungapped_sequence = sequence_to_search[start_index:end_index]

                        if is_rc:
                            pos_in_original = total_len - end_index
                        else:
                            pos_in_original = start_index

                        if ungapped_to_gapped is None:
                            ungapped_to_gapped = np.flatnonzero(seq_record.seq_np != 45)

                        if pos_in_original + match_len - 1 >= len(ungapped_to_gapped):
                            continue
                        gapped_start = int(ungapped_to_gapped[pos_in_original])
                        gapped_end = int(ungapped_to_gapped[pos_in_original + match_len - 1])

                        gapped_sequence = seq_gapped[gapped_start:gapped_end + 1]

                        vectors, conservation_value, aligned_sequences, genomic_start, genomic_end, ref_chromosome = compute_vectors_and_conservation(
                            block,
                            ref_seq_record,
                            genome_ids,
                            gapped_start,
                            gapped_end + 1,
                            ref_genome_id,
                            local_genome_names,
                            pos_in_original,
                            match_len
                        )

                        motif_hit_data = {
                            "hit_id": f"regex_{genome_id}_{regex_identifier}{hit_id_suffix}_{pos_in_original}",
                            "reference_genome": ref_genome_id,
                            "chromosome": ref_chromosome,
                            "genomic_start": genomic_start,
                            "genomic_end": genomic_end,
                            "strand": sequence_strand,
                            "motif_type": "regex",
                            "motif_name": regex_identifier,
                            "motif_length": match_len,
                            "gapped_motif": gapped_sequence,
                            "gapped_start": gapped_start,
                            "gapped_end": gapped_end,
                            "ungapped_start": pos_in_original,
                            "ungapped_end": pos_in_original + match_len - 1,
                            "motif": ungapped_sequence,
                            "motif_strand": motif_strand,
                            "conservation": {
                                "value": float(round(conservation_value, 2)),
                                "other_genomes": aligned_sequences
                            }
                        }

                        key = (genome_id, regex_identifier)
                        if key not in genome_files:
                            genome_files[key] = []
                        
                        try:
                            genome_files[key].append(fast_dumps(motif_hit_data))
                            if len(genome_files[key]) >= 10000:
                                genome_dir = os.path.join('tmp', genome_id)
                                output_file = os.path.join(
                                    genome_dir,
                                    f'{output_identifier}_{regex_identifier}_motif_hits_process_{process_id}.tmp'
                                )
                                io_executor.submit(write_hits_to_disk, genome_dir, output_file, list(genome_files[key]))
                                genome_files[key].clear()
                        except Exception as e:
                            logger.error(f"Error buffering motif hit: {e}", exc_info=True)

            else:
                logger.error(f"Unknown search type: {search_type}")
                sys.exit(1)

    except Exception as e:
        logger.error(f"Error processing block {block_no}: {e}", exc_info=True)
        sys.exit(1)


def _initialize_automaton(patterns):
    """Initialize an Aho-Corasick automaton for k-mer searching."""
    A = ahocorasick.Automaton()
    for idx, kmer_data in enumerate(patterns):
        kmer = kmer_data['kmer']
        kmer_data['index'] = idx
        A.add_word(kmer, (idx, kmer_data))
    A.make_automaton()
    return A


def process_file_chunk(maf_file,
                       start_pos,
                       end_pos,
                       search_type,
                       patterns,
                       genome_ids,
                       search_in,
                       reverse_complement,
                       process_id,
                       unique_genome_names,
                       output_identifier):
    """Process a chunk of the MAF file and write results to .tmp files."""
    genome_files = {}
    block_no = 0
    local_genome_names = set()
    A = _initialize_automaton(patterns) if search_type == "kmer" else None
    
    # Pre-allocate a large pool of SimpleSeqRecord objects to avoid 
    # instantiating thousands of objects during file iteration loops.
    record_pool = [SimpleSeqRecord() for _ in range(500)]
    io_executor = ThreadPoolExecutor(max_workers=2)

    try:
        with open(maf_file, 'rb') as handle:
            handle.seek(start_pos)
            chunk_data = handle.read(end_pos - start_pos).decode('utf-8')

        # Split blocks exactly by '\na' indicating a new MAF block
        blocks = chunk_data.split('\na')
        
        for idx, block_data in enumerate(blocks):
            if not block_data.strip():
                continue
                
            # split('\na') removes the 'a', so restore it for all blocks EXCEPT the first
            # (unless the first block actually had a newline before it, which is rare at boundaries)
            block_string = 'a' + block_data if idx > 0 else block_data

            block = parse_block_from_string(block_string, record_pool)
            if block:
                search_patterns_in_block(
                    block,
                    search_type,
                    A,
                    patterns,
                    genome_ids,
                    block_no,
                    search_in,
                    genome_files,
                    reverse_complement,
                    process_id,
                    local_genome_names,
                    output_identifier,
                    io_executor
                )
            block_no += 1

    except Exception as e:
        logger.error(f"Error in process {process_id}: {e}", exc_info=True)
    finally:
        # Flush any remaining built-up motif arrays into their respective tmp files
        for key, buffer in genome_files.items():
            if buffer:
                # Based on the key, rebuild the identifier
                # For k-mer: key is genome_id
                # For regex/pwm: key is (genome_id, motif_identifier)
                if isinstance(key, tuple):
                    genome_id, motif_identifier = key
                    filename = f'{output_identifier}_{motif_identifier}_motif_hits_process_{process_id}.tmp'
                else:
                    genome_id = key
                    filename = f'{output_identifier}_motif_hits_process_{process_id}.tmp'
                
                genome_dir = os.path.join('tmp', genome_id)
                os.makedirs(genome_dir, exist_ok=True)
                output_file = os.path.join(genome_dir, filename)
                
                try:
                    io_executor.submit(write_hits_to_disk, genome_dir, output_file, list(buffer))
                except Exception as e:
                    logger.error(f"Error flushing genome file: {e}", exc_info=True)

        io_executor.shutdown()
        unique_genome_names.extend(local_genome_names)


def generate_bed(tmp_dir, bed_filename):
    """
    Parse all .tmp files in tmp_dir, collect motif hits,
    and write them to a single BED (TSV) file with columns:
    chrom, chromStart, chromEnd, name, score, strand

    Coordinates remain 0-based inclusive.
    """
    bed_entries = []
    # Recursively find all .tmp files
    for root, dirs, files in os.walk(tmp_dir):
        for f in files:
            if f.endswith('.tmp'):
                tmp_file_path = os.path.join(root, f)
                try:
                    with open(tmp_file_path, 'r') as fh:
                        for line in fh:
                            line = line.strip()
                            if not line:
                                continue
                            try:
                                data = fast_json.loads(line)
                                chrom = data.get('chromosome', 'unknown')
                                chrom_start = data.get('genomic_start', 0)  # 0-based inclusive
                                chrom_end = data.get('genomic_end', 0)+1    # 0-based non inclusive (open)
                                name = data.get('motif_name', 'motif')
                                score_val = data.get('conservation', {}).get('value', 0.0)
                                if isinstance(score_val, (float, int)):
                                    score_val = round(score_val, 2)
                                else:
                                    score_val = 0.0
                                strand = data.get('strand', '+')

                                # Write to BED file
                                bed_entries.append(
                                    f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score_val}\t{strand}"
                                )
                            except ValueError:
                                logger.error(f"Invalid JSON in {tmp_file_path}")
                except Exception as e:
                    logger.error(f"Error reading .tmp file {tmp_file_path}: {e}")

    if bed_entries:
        with open(bed_filename, 'w') as bed_out:
            for entry in bed_entries:
                bed_out.write(entry + "\n")


def merge_results(tmp_dir, output_identifier, detailed_report):
    """
    Merge .tmp files into per-genome/per-motif merged JSON if --detailed_report is used.
    If not detailed_report, skip the JSON merging.
    Return a dict { (genome_id, motif_identifier) : list_of_hits } for CSV usage if needed.
    """
    if not detailed_report:
        # If user didn't request detailed report, skip merging JSON
        return {}

    logger.info("Merging .tmp results into JSON (detailed report).")
    genome_dirs = [d for d in os.listdir(tmp_dir) if os.path.isdir(os.path.join(tmp_dir, d))]
    all_merged_data = {}  # (genome_id, motif_id) -> list of hits

    for genome_id in genome_dirs:
        genome_path = os.path.join(tmp_dir, genome_id)
        tmp_files = [x for x in os.listdir(genome_path) if x.endswith('.tmp')]

        files_by_motif = {}
        for f in tmp_files:
            # Example:  <output_identifier>_<motif_identifier>_motif_hits_process_<pid>.tmp
            match = re.match(rf'{output_identifier}_(.+?)_motif_hits_process_\d+\.tmp', f)
            if match:
                motif_id = match.group(1)
                files_by_motif.setdefault(motif_id, []).append(f)
            else:
                # Could be kmer with no motif_id in the filename
                if f.startswith(f'{output_identifier}_motif_hits_process_'):
                    motif_id = output_identifier
                    files_by_motif.setdefault(motif_id, []).append(f)

        for motif_identifier, file_list in files_by_motif.items():
            merged_data = []
            # Sort by process ID if present
            file_list.sort(
                key=lambda x: int(re.findall(r'\d+', x)[-1]) if re.findall(r'\d+', x) else 0
            )
            for tmp_file in file_list:
                file_path = os.path.join(genome_path, tmp_file)
                if os.path.exists(file_path):
                    with open(file_path, 'r') as fh:
                        for line in fh:
                            line = line.strip()
                            if not line:
                                continue
                            try:
                                data = fast_json.loads(line)
                                merged_data.append(data)
                            except ValueError:
                                logger.error(f"Error decoding JSON from {file_path}")

            if merged_data:
                all_merged_data.setdefault((genome_id, motif_identifier), []).extend(merged_data)

    # Write them out to JSON (in current directory, not tmp) if not empty
    for (genome_id, motif_identifier), data_list in all_merged_data.items():
        if not data_list:
            continue
        out_file_name = f"{genome_id}_{output_identifier}_{motif_identifier}_motif_hits.json"
        try:
            with open(out_file_name, 'w') as f:
                f.write(fast_dumps(data_list))
            logger.info(f"Wrote merged JSON: {out_file_name}")
        except Exception as e:
            logger.error(f"Error writing merged results {out_file_name}: {e}")

    return all_merged_data


def generate_csv(all_merged_data, unique_genome_names, output_identifier):
    """
    Generate CSV files if --detailed_report is used.
    We have the dictionary from merge_results: {(genome_id, motif_identifier): [hits]}
    Write them in current directory as genomeID_outputIdentifier_motifID_motif_hits.csv

    Coordinates in the CSV are also 0-based inclusive.
    """
    if not all_merged_data:
        logger.info("No merged data found; skipping CSV generation.")
        return

    logger.info("Generating CSV files (detailed report).")
    unique_genome_names = set(unique_genome_names)

    for (genome_id, motif_identifier), data_list in all_merged_data.items():
        # Build a CSV
        out_file_name = f"{genome_id}_{output_identifier}_{motif_identifier}_motif_hits.csv"
        all_motif_hits = {}

        for motif_hit in data_list:
            hit_id = motif_hit['hit_id']
            all_motif_hits[hit_id] = motif_hit

        fieldnames = ['motif_hit_info', 'motif_name', 'motif_type', 'motif_strand', 'score'] + list(unique_genome_names)
        try:
            with open(out_file_name, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

                for hit_id, motif_hit in all_motif_hits.items():
                    motif_hit_info = (
                        f"{motif_hit['chromosome']}:{motif_hit['genomic_start']}-{motif_hit['genomic_end']},"
                        f"{motif_hit['motif']}"
                    )
                    row = {
                        'motif_hit_info': motif_hit_info,
                        'motif_name': motif_hit['motif_name'],
                        'motif_type': motif_hit['motif_type'],
                        'motif_strand': motif_hit['motif_strand'],
                        'score': round(motif_hit.get('conservation', {}).get('value', 0.0), 2),
                    }
                    for gname in unique_genome_names:
                        row[gname] = None

                    other_genomes = motif_hit.get('conservation', {}).get('other_genomes', [])
                    for genome_data in other_genomes:
                        gname = genome_data['genome_id']
                        if gname in row:
                            vector = genome_data.get('vector')
                            conservation_pct = genome_data.get('conservation')
                            row[gname] = f"{vector},{conservation_pct}"

                    writer.writerow(row)
            logger.info(f"CSV file generated: {out_file_name}")
        except Exception as e:
            logger.error(f"Error writing CSV file {out_file_name}: {e}")


def purge_directory(directory):
    if os.path.exists(directory):
        logger.info(f"Purging {directory}...")
        shutil.rmtree(directory)
        logger.info(f"{directory} has been purged.")
    else:
        logger.info(f"{directory} does not exist.")


def main():
    parser = argparse.ArgumentParser(
        description="Search uncompressed MAF files for specific motifs using regex patterns, K-mers, or PWMs."
    )

    parser.add_argument('maf_file', help="Path to the uncompressed MAF file.")
    parser.add_argument('--genome_ids', help="Path to the genome IDs file (optional).")
    parser.add_argument('--search_in', default='reference', help="Genome in which to search motifs (default: reference).")
    parser.add_argument('--reverse_complement', choices=['yes', 'no'], default='no',
                        help="Search motifs on both strands if 'yes'. (Default: no)")
    parser.add_argument('--pvalue_threshold', type=float, default=1e-4,
                        help="P-value threshold for PWM matches. (Default: 1e-4)")
    parser.add_argument('--processes', type=int, default=1,
                        help="Number of parallel processes to use. (Default: 1)")
    parser.add_argument('--background_frequencies', nargs=4, type=float,
                        metavar=('A_FREQ', 'C_FREQ', 'G_FREQ', 'T_FREQ'),
                        help="Background nucleotide frequencies for A, C, G, T. Must sum to 1.")
    parser.add_argument('--purge_results_dir', action='store_true',
                        help="Purge the 'tmp' directory before running.")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help="Enable verbose logging to debug.log.")
    parser.add_argument('--detailed_report', '-d', action='store_true',
                        help="Generate merged JSON and CSV output in current directory.")

    # Mutually exclusive group for the search type
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--regexes', nargs='+', help="Regex patterns to search.")
    group.add_argument('--kmers', help="Path to a file of K-mers (one per line).")
    group.add_argument('--jaspar_file', help="Path to a JASPAR-format PWM file.")

    args = parser.parse_args()

    # Configure logging based on --verbose
    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler('debug.log', mode='w')]
        )
    else:
        # Minimal logging
        logging.basicConfig(level=logging.CRITICAL)

    logger.info("Starting processing...")

    maf_file = args.maf_file
    genome_ids_file = args.genome_ids
    search_in = args.search_in
    reverse_complement = (args.reverse_complement == 'yes')
    pvalue_threshold = args.pvalue_threshold
    num_processes = args.processes

    tmp_dir = os.path.join(os.getcwd(), 'tmp')
    if args.purge_results_dir:
        purge_directory(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)

    genome_ids = load_genome_ids(genome_ids_file)

    # Process background frequencies
    background_frequencies = None
    if args.background_frequencies:
        total = sum(args.background_frequencies)
        if not np.isclose(total, 1.0):
            logger.error("The background frequencies must sum to 1.")
            print("Failure!")
            sys.exit(1)
        nucleotides = ['A', 'C', 'G', 'T']
        background_frequencies = dict(zip(nucleotides, args.background_frequencies))
        logger.info(f"Using background frequencies: {background_frequencies}")

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
                regex_patterns.append({
                    'pattern': pattern,
                    'regex': regex,
                    'identifier': regex_identifier,
                    'is_reverse_complement': False
                })
                if reverse_complement:
                    regex_patterns.append({
                        'pattern': pattern,
                        'regex': regex,
                        'identifier': regex_identifier,
                        'is_reverse_complement': True
                    })
            except re.error as e:
                logger.error(f"Invalid regex pattern '{regex}': {e}")
                print("Failure!")
                sys.exit(1)
        patterns = regex_patterns
        logger.info(f"Loaded {len(patterns)} regex patterns.")
    elif args.kmers:
        search_type = 'kmer'
        kmer_patterns = []
        kmers_file_name = os.path.splitext(os.path.basename(args.kmers))[0]
        output_identifier = kmers_file_name
        try:
            with open(args.kmers, 'r') as f:
                kmers_list = list(set([line.strip().upper() for line in f if line.strip()]))
            for kmer in kmers_list:
                kmer_patterns.append({
                    'kmer': kmer,
                    'identifier': kmer,
                    'original_kmer': kmer,
                    'is_reverse_complement': False
                })
                if reverse_complement:
                    rc_kmer = str(Seq(kmer).reverse_complement())
                    kmer_patterns.append({
                        'kmer': rc_kmer,
                        'identifier': kmer,
                        'original_kmer': kmer,
                        'is_reverse_complement': True
                    })
            # Remove duplicates
            kmer_patterns = list({(d['kmer'], d['is_reverse_complement']): d for d in kmer_patterns}.values())
            patterns = kmer_patterns
            logger.info(f"Loaded {len(patterns)} K-mers (including reverse complements).")
        except Exception as e:
            logger.error(f"Error loading K-mers from {args.kmers}: {e}", exc_info=True)
            print("Failure!")
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

                    if background_frequencies:
                        pssm.background = background_frequencies
                    elif motif.background:
                        pssm.background = motif.background
                    else:
                        pssm.background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

                    threshold_score = compute_threshold_score(pssm, pvalue_threshold)
                    pwm_motifs.append({
                        'pwm': pssm,
                        'name': motif.name,
                        'threshold': threshold_score,
                        'identifier': motif_identifier,
                        'is_reverse_complement': False
                    })

                    if reverse_complement:
                        pssm_rc = pssm.reverse_complement()
                        pssm_rc.background = pssm.background
                        threshold_rc = compute_threshold_score(pssm_rc, pvalue_threshold)
                        pwm_motifs.append({
                            'pwm': pssm_rc,
                            'name': motif.name,
                            'threshold': threshold_rc,
                            'identifier': motif_identifier,
                            'is_reverse_complement': True
                        })

            patterns = pwm_motifs
            logger.info(f"Loaded {len(patterns)} PWMs (including reverse complements).")
        except Exception as e:
            logger.error(f"Error loading PWMs from {args.jaspar_file}: {e}", exc_info=True)
            print("Failure!")
            sys.exit(1)
    else:
        logger.error("One of --regexes, --kmers, or --jaspar_file must be provided.")
        print("Failure!")
        sys.exit(1)

    chunk_ranges = assign_file_chunks(maf_file, num_processes)

    manager = multiprocessing.Manager()
    unique_genome_names = manager.list()
    processes = []

    for i, (start_pos, end_pos) in enumerate(chunk_ranges):
        p = multiprocessing.Process(
            target=process_file_chunk,
            args=(
                maf_file, start_pos, end_pos, search_type, patterns,
                genome_ids, search_in, reverse_complement, i,
                unique_genome_names, output_identifier
            )
        )
        processes.append(p)
        p.start()
        logger.info(f"Started process {i} with PID {p.pid}")

    for i, p in enumerate(processes):
        p.join()
        logger.info(f"Process {i} with PID {p.pid} has completed.")

    unique_genome_names = list(set(unique_genome_names))
    logger.info(f"Unique genome names encountered: {unique_genome_names}")

    # Always generate a BED file from the .tmp outputs
    bed_filename = f"{output_identifier}_motif_hits.bed"
    generate_bed(tmp_dir, bed_filename)
    logger.info(f"BED file generated: {bed_filename}")

    # If --detailed_report, also merge to JSON & CSV in current dir
    all_merged_data = merge_results(tmp_dir, output_identifier, args.detailed_report)
    if args.detailed_report:
        generate_csv(all_merged_data, unique_genome_names, output_identifier)

    # Purge tmp directory
    purge_directory(tmp_dir)

    logger.info("Completed processing.")
    print("Success!")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("Failure!")
        sys.exit(1)
