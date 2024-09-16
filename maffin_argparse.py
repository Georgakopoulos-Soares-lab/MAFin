import sys

import gzip

import csv

import re

import argparse

from Bio import AlignIO

import pandas as pd

import matplotlib.pyplot as plt



def load_genome_ids(genome_ids_file):

    """Load genome IDs from a given file."""

    if genome_ids_file:

        print(f"Loading genome IDs from {genome_ids_file}", flush=True)

        with open(genome_ids_file, 'r') as file:

            return file.read().splitlines()

    return None



def find_genome_ids_in_maf(maf_file):

    """Find all unique genome IDs present in the MAF file."""

    genome_ids_found = set()

    with gzip.open(maf_file, "rt") as handle:

        alignment_blocks = AlignIO.parse(handle, "maf")

        for block in alignment_blocks:

            for seq in block:

                genome_id = seq.id.split('.')[0]

                genome_ids_found.add(genome_id)

    return list(genome_ids_found)



def compute_similarity_vector(kmer, block, ref_start_index, genome_ids, extend):

    """Compute similarity vector for a kmer within an alignment block,

    comparing other sequences to the specific reference segment where the kmer was found,

    with gaps removed from all sequences."""

    kmer_length = len(kmer)

    extended_length = kmer_length + 20 if extend else kmer_length

    

    # Initialize vectors and conservation percentages with "NA"

    vectors = {genome_id: 'NA' for genome_id in genome_ids}

    conservation_percentages = {genome_id: "NA" for genome_id in genome_ids}



    # Extract the reference sequence and remove gaps

    ref_seq = str(block[0].seq).replace('-', '')

    

    # Calculate the start and end indices for the reference segment

    start_idx = max(0, ref_start_index - 10) if extend else ref_start_index

    end_idx = start_idx + extended_length



    # Ensure the end index does not exceed the reference sequence length

    ref_segment = ref_seq[start_idx:end_idx]



    # Handle the reference genome

    ref_genome_id = block[0].id.split('.')[0]

    vectors[ref_genome_id] = '1' * extended_length

    conservation_percentages[ref_genome_id] = "100.00%"



    for seq in block[1:]:

        trimmed_id = seq.id.split('.')[0]

        if trimmed_id in vectors:

            ungapped_seq = str(seq.seq).replace('-', '')

            seq_segment = ungapped_seq[start_idx:end_idx]

            vector = ['0'] * extended_length

            match_count = 0



            # Compare each position

            for i in range(min(len(ref_segment), len(seq_segment))):

                if seq_segment[i] == ref_segment[i]:

                    vector[i] = '1'

                    if start_idx + i >= ref_start_index and start_idx + i < ref_start_index + kmer_length:

                        match_count += 1



            vectors[trimmed_id] = ''.join(vector)

            percentage = (match_count / kmer_length) * 100 if kmer_length > 0 else 0

            conservation_percentages[trimmed_id] = f"{percentage:.2f}%"



    return vectors, conservation_percentages



def calculate_average_conservation(kmer_conservation):

    """Calculate the average conservation percentage for a kmer, ignoring 'NA' values."""

    total_percentage = 0

    count = 0

    for percentage in kmer_conservation.values():

        if percentage != 'NA':

            total_percentage += float(percentage.replace('%', ''))

            count += 1

    return total_percentage / count if count > 0 else 0



def search_kmers_in_block(block, regex_patterns, genome_ids, csv_writers, regex_avg_conservation):

    """Search for kmers matching regex patterns in the reference sequence of the alignment block."""

    reference = block[0]

    sequence = str(reference.seq).replace('-', '')

    seq_start = reference.annotations['start']



    for regex_pattern in regex_patterns:

        for match in re.finditer(regex_pattern, sequence):

            start = match.start()

            end = match.end()

            kmer = match.group()

            ref_index = start

            extend = (ref_index - 10 >= 0) and (ref_index + len(kmer) + 10 <= len(sequence))

            vectors, conservation_percentages = compute_similarity_vector(kmer, block, ref_index, genome_ids, extend)



            # Calculate average conservation for this kmer and store

            avg_conservation = calculate_average_conservation(conservation_percentages)

            regex_avg_conservation[regex_pattern].append(avg_conservation)



            line = f"{reference.id.split('.')[1]},{seq_start + start},{seq_start + end},{kmer}"

            row = [line] + [f"{vectors[genome]} ({conservation_percentages[genome]})" for genome in genome_ids] + ["YES" if extend else "NO"]



            csv_writers[regex_pattern].writerow(row)



def main(maf_file, regex_patterns, genome_ids_file=None):

    print("Starting..", flush=True)



    # Find all genome IDs in the MAF file

    all_genome_ids = find_genome_ids_in_maf(maf_file)

    print(f"Found {len(all_genome_ids)} unique genome IDs in the MAF file.", flush=True)

    

    # Load the genome IDs to consider from genome_ids.txt, if provided

    genome_ids_to_use = load_genome_ids(genome_ids_file)

    if genome_ids_to_use:

        # Keep only the genome IDs that are present in the genome_ids.txt file

        genome_ids = [genome for genome in all_genome_ids if genome in genome_ids_to_use]

        print(f"Considering {len(genome_ids)} genome IDs based on the genome_ids.txt file.", flush=True)

    else:

        genome_ids = all_genome_ids

        print(f"Considering all {len(genome_ids)} genome IDs found in the MAF file.", flush=True)



    # Create output CSV files for each regex pattern

    csv_writers = {}

    regex_avg_conservation = {regex: [] for regex in regex_patterns}

    for regex_pattern in regex_patterns:

        safe_regex_name = re.sub(r'\W+', '_', regex_pattern)

        output_file = f"output_{safe_regex_name}.csv"

        out_file = open(output_file, mode='w', newline='')

        csv_writer = csv.writer(out_file)

        # Write the header with all genome IDs

        header = ["Input Line"] + [f"{genome_id} (Vector, %Conservation)" for genome_id in genome_ids] + ["+/- 10 BP"]

        csv_writer.writerow(header)

        csv_writers[regex_pattern] = csv_writer



    # Process the MAF file again to search for kmers and write to CSV

    with gzip.open(maf_file, "rt") as handle:

        print(f"Opened MAF file: {maf_file}", flush=True)

        alignment_blocks = AlignIO.parse(handle, "maf")

        for block in alignment_blocks:

            search_kmers_in_block(block, regex_patterns, genome_ids, csv_writers, regex_avg_conservation)



    # Calculate overall average conservation for each regex

    overall_avg_conservation = {}

    for regex_pattern, kmer_conservations in regex_avg_conservation.items():

        if kmer_conservations:

            overall_avg_conservation[regex_pattern] = sum(kmer_conservations) / len(kmer_conservations)

        else:

            overall_avg_conservation[regex_pattern] = 0



    # Plot the average conservation for each regex

    plt.bar(overall_avg_conservation.keys(), overall_avg_conservation.values())

    plt.xlabel("Regex Patterns")

    plt.ylabel("Average Conservation (%)")

    plt.title("Average Conservation for Each Regex Pattern")

    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()

    plt.savefig('average_conservation_plot.png')  # Save as PNG

    plt.savefig('average_conservation_plot.pdf')  # Optionally, save as PDF

    plt.show()



    print("Completed processing.", flush=True)



if __name__ == "__main__":

    # Create an argument parser

    parser = argparse.ArgumentParser(description="Search MAF files for specific kmers using regex patterns.")

    

    # Required arguments

    parser.add_argument('--maf_file', required=True, help="Path to the MAF file.")

    parser.add_argument('--regexes', required=True, nargs='+', help="Regex patterns to search for in the MAF file.")

    

    # Optional arguments

    parser.add_argument('--genome_ids', help="Path to the genome IDs file (optional). If not provided, all genomes are used.")



    # Parse the arguments

    args = parser.parse_args()



    # Validate that MAF file and regexes are provided

    if not args.maf_file:

        print("Error: --maf_file must be provided.")

        sys.exit(1)

    

    if not args.regexes:

        print("Error: --regexes must be provided.")

        sys.exit(1)



    # Run the main function with the parsed arguments

    main(args.maf_file, args.regexes, args.genome_ids)


