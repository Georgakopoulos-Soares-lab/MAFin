<!-- MAFin Logo -->
![Alt text](https://i.postimg.cc/V5QnxMWp/logo-no-background.png "")

# MAFin: Motif Detection in Multiple Alignment Files

## Description

This CLI tool allows you to search MAF (Multiple Alignment Format) files for specific motifs using regular expression patterns, K-mers, or Position Weight Matrices (PWMs). It supports searching within reference genomes or across all genomes and can handle both strands if specified. The tool is optimized for performance with support for parallel processing.

## Features

- **Motif Search Options**: Search using regex patterns, K-mers, or PWMs (JASPAR format).
- **Flexible Genome Selection**: Choose to search within the reference genome or across all genomes.
- **Strand Orientation**: Optionally search both strands by considering reverse complements.
- **Parallel Processing**: Utilize multiple processes to speed up the search.
- **Output Formats**: Generate detailed JSON and CSV reports of the motifs found.

## Installation

### Prerequisites

- Python 3.6 or higher
- `pip` package manager
-  Biopython, pyahocorasick ( installed via requirements.txt )

### Using Virtual Environment

1. **Clone the Repository**

   ```bash
   git clone https://github.com/Georgakopoulos-Soares-lab/MAFin
   cd MAFin
    ```

2. **Install package inside conda environment**
   ```bash
   conda create -n mafin_env python=3.8
   conda init && conda activate mafin_env
   pip install -e .
   ```
   or through PyPI without cloning
   ```bash
   pip install MAFin   
   ```

## Usage
   ```bash
   MAFin --maf_file path/to/file.maf \
                           --regexes "CTGCCCGCA" "AGT" \
                           --search_in reference \
                           --reverse_complement no \
                           --pvalue_threshold 1e-4 \
                           --processes 4
   ```


### Command-Line Arguments

    --maf_file: (Required) Path to the uncompressed MAF file.
    --genome_ids: Path to the genome IDs file (optional). If not provided, all genomes are used.
    --search_in: Specify whether to search motifs in the 'reference' genome or 'all' genomes. Default is 'reference'.
    --reverse_complement: Specify whether to search motifs on both strands. Choices are 'yes' or 'no'. Default is 'no'.
    --pvalue_threshold: P-value threshold for PWM matches. Default is 1e-4.
    --processes: Number of parallel processes to use. Default is 1.
    --background_frequencies: Background nucleotide frequencies for A, C, G, T. They must sum to 1. ( defaults 0.25 0.25 0.25 0.25)
    --purge_results_dir: Purge the results directory before running the script.
    
### Search Types (Mutually Exclusive)

    --regexes: Regex patterns to search for in the MAF file.
    --kmers: Path to the K-mers file (one per line).
    --jaspar_file: Path to the PWM file in JASPAR format.

<b>Note:</b> You must specify exactly one of the search types.

## Example Usage


### Searching with Regex Patterns
```bash
python maf_motif_search.py --maf_file data/sample.maf \
                           --regexes "CTG[CC]+CGCA" "AGT" \
                           --search_in all \
                           --reverse_complement yes \
                           --pvalue_threshold 0.001 \
                           --processes 8 \
                           --purge_results_dir
```

### Searching with K-mers
```bash
                  MAFin    --maf_file data/sample.maf \
                           --kmers data/kmers.txt \
                           --search_in reference \
                           --pvalue_threshold 1e-5 \
                            --background_frequencies 0.25 0.15 0.35 0.25 \
                           --processes 2
```

### Searching with PWM (JASPAR format)
```bash
                     MAFin --maf_file data/sample.maf \
                           --jaspar_file data/motif.jaspar \
                           --search_in all \
                           --processes 4
```

## Example Output

### JSON Output

After running the tool, a JSON file is generated containing detailed information about the motifs found.
```json
{
    "num_blocks": 16393,
    "count_of_genomes": 6,
    "genome_ids": [
        "HG00621#2",
        "HG02723#2",
        "HG01175#2",
        "NA20129#1",
        "HG03098#1",
        "HG01123#1"
    ],
    "motifs_found": 89,
    "motifs": [
        {
            "hit_id": "pssm_THAP1_1",
            "block_no": 433,
            "motif_type": "pssm",
            "motif_length": 9,
            "motif": "CTGCCCGCA",
            "gapped_start": 1038,
            "gapped_end": 1046,
            "ungapped_start": 1038,
            "ungapped_end": 1046,
            "gapped_sequence": "CTGCCCGCA",
            "ungapped_sequence": "CTGCCCGCA",
            "conservation": {
                "value": 100.0,
                "other_genomes": [
                    {
                        "genome_id": "CHM13",
                        "aligned_sequence_gapped": "CTGCCCGCA",
                        "vector": "111111111",
                        "conservation": "100.00%",
                        "gapped_start": 1038,
                        "gapped_end": 1046
                    },
                    {
                        "genome_id": "GRCh38",
                        "aligned_sequence_gapped": "CTGCCCGCA",
                        "vector": "111111111",
                        "conservation": "100.00%",
                        "gapped_start": 1038,
                        "gapped_end": 1046
                    }
                    // Additional genome entries...
                ]
            }
        }
        // Additional motif entries...
    ]
}

```

### CSV Output

The csv file generated contains the coordinates of the motif found along with its name , type , strand direction and conservation. Each motif_hit is then followed by a number of columns for each aligned sequence in the MAF block where the hit was found. The cells are then filled with an entry, <b>similarity_vector, conservation_score</b>
<br>
e.g. <b>11110110000,54.55%</b>
 	

| motif_hit_info            | motif_name | motif_type | motif_strand | score               | fake_species_48             | fake_species_45             | fake_species_6              | fake_species_54             | fake_species_65             |
|---------------------------|------------|------------|--------------|---------------------|-----------------------------|-----------------------------|-----------------------------|-----------------------------|-----------------------------|
| chr16:9722-9732,ACCTGTGGTCT | MA0002.2   | pwm        | +            | 10.025236129760742   | 11110110000,54.55%          | 11111111111,100.00%         | 11111101111,90.91%          | 11011111111,90.91%          | 10110111011,72.73%          |
| chr46:33765-33775,TAACCACAGAA | MA0002.2   | pwm        | +            | 11.518963813781738   | 11111111111,100.00%         |                             | 11111111011,90.91%          | 11111111111,100.00%         | 11011011110,72.73%          |
| chr69:51945-51955,TGTTGTGGTTA | MA0002.2   | pwm        | +            | 10.343249320983887   | 11111111111,100.00%         |                             | 01111111111,90.91%          | 11111111111,100.00%         | 11111011110,81.82%          |
| chr80:61982-61992,CAACCACAGGC | MA0002.2   | pwm        | +            | 10.74906063079834    |                             | 10111111111,90.91%          | 10111111011,81.82%          |                             | 11111111101,90.91%          |


## Workflow Diagram
[![maffin-drawio-1.png](https://i.postimg.cc/nLsbvLmz/maffin-drawio-1.png)](https://postimg.cc/Mc8Fwqd2)

## Conservation Logic in MAFin

MAFin is a tool that helps identify conserved motifs across multiple genomes by comparing aligned sequences in MAF files. These motifs can be found using one of three methods:

- **Regular Expression (regex)**: A pattern search within the sequence.
- **Position Weight Matrix (PWM)**: A motif search using a predefined PWM.
- **K-mer Input**: Using specific k-mers to search the sequences.

In the following examples, we will assume the motif being searched for is `ATCG`.

### Conservation Process

Once the motif has been found, MAFin calculates the conservation across aligned sequences. The process is based on comparing the **ungapped sequences** within the alignment, while still respecting the gaps in the MAF file to maintain the alignment structure. The comparison results in a **similarity vector** that matches the true length of the motif, excluding positions where both the reference and compared sequences have gaps.

The conservation logic also includes calculating a **conservation percentage**, which tells us how conserved the motif is across genomes, based on the number of matches relative to the length of the motif.

### Conservation Logic Rules

1. **Gap in Both Sequences (Skip)**:
   - If both the reference genome and the compared genome have a gap at the same position, we skip this position.

2. **Match (1)**:
   - If both genomes have the same base pair at the position, it's considered a match, and we mark it as `1`.

3. **Mismatch (0)**:
   - If one genome has a gap and the other has a base, or if the bases are different, it is considered a mismatch, and we mark it as `0`.

### Conservation Score Calculation

The conservation percentage is calculated based on the number of matches relative to the length of the motif (excluding positions where both sequences have gaps). The formula for conservation is:

<!-- You might want to add the actual formula here if it's missing -->

### Simple Case Example: Matching Sequences

Let’s start with a simple example where the motif `ATCG` is found in both the reference genome and the compared genome with identical alignment and gaps.

**Motif**: `ATCG`

**Reference Genome**: `A - T C G`  
**Compared Genome**:  `A - T C G`

In this case, the sequences are exactly the same, with gaps aligned in the same positions. MAFin will compare the ungapped bases to maintain the alignment structure.

**Step-by-Step Comparison**:

| Position | Ref Base | Compared Base | Result | Vector |
|----------|----------|---------------|--------|--------|
| 1        | A        | A             | Match  | 1      |
| 2        | -        | -             | Skip   |        |
| 3        | T        | T             | Match  | 1      |
| 4        | C        | C             | Match  | 1      |
| 5        | G        | G             | Match  | 1      |

**Similarity Vector**:  
Since we skip the gap in position 2, the resulting similarity vector has a length of 4, matching the true length of the motif (`ATCG`):

<!-- You might want to include the actual similarity vector here if it's missing -->

**Conservation Percentage**:  
The total positions compared (excluding gaps) is 4, and all of them are matches. Therefore, the conservation score is: 100%

<!-- You might want to include the actual conservation percentage calculation here if it's missing -->

### Genomic Coordinates for Matches

Along with the similarity vector, MAFin provides the genomic coordinates of the motif. If the motif starts at position 1000 in the reference genome and spans 4 ungapped bases (`A T C G`), the start and end positions would be: start: 1000, end:1003

MAFin provides these coordinates for the motif in both the reference genome and the aligned genomes.

### Example 2: Mismatched Sequences with Gaps

Now, let’s look at a more complex case where gaps and mismatches occur between the reference and compared sequences.

**Motif**: `ATCG`

**Reference Genome**: `A - T C G`  
**Compared Genome**: `A T - C G`

Here, the sequences differ, with gaps in different positions. MAFin will again compare the ungapped bases to produce a similarity vector of length 4 (the length of the motif).

**Step-by-Step Comparison**:

| Position | Ref Base | Compared Base | Result   | Vector |
|----------|----------|---------------|----------|--------|
| 1        | A        | A             | Match    | 1      |
| 2        | -        | T             | Mismatch | 0      |
| 3        | T        | -             | Mismatch | 0      |
| 4        | C        | C             | Match    | 1      |
| 5        | G        | G             | Match    | 1      |

**Similarity Vector**:  
Again, we skip the gap in position 2 of the reference sequence, and the resulting similarity vector still has a length of 4: [1,0,0,1]. Note here that the last element is missing.


**Conservation Percentage**:  
The total positions compared (excluding gaps) is 4, and 2 out of 4 positions are matches. Therefore, the conservation score is: 50%

<!-- You might want to include the actual conservation percentage calculation here if it's missing -->

### Genomic Coordinates for Matches and Mismatches

As with the first example, MAFin provides genomic coordinates for the motif. If the motif starts at position 1000 in the reference genome and spans 4 ungapped positions, the coordinates are: start:1000, end:1003

These coordinates, along with the similarity vector and conservation score, allow users to easily trace the conserved motifs across genomes.

---

## Reverse Complement Example

Searcing for a motif on the reverse strand generally involves a search for the reverse complement of a sequence.

### Searching reverse strand through Regex Patterns

**Reference Sequence**: `ATCGGCA`
**Regular Expression**: `C{2}G` ( Two times C followed by G )

Due to the complex nature of Regex patterns reversing and complementing the expression is impossible. Thus we reverse and complement the reference sequence and then search for the pattern there.

Reverse complement match of CCG in sequence: TG**CCG**AT , results in genomic coordinates 3,5:
| 1 | 2 | 3 | 4   | 5 |6|7|
|----------|----------|---------------|----------|--------|-|-|
| T        | G        | **C**             | **C**    | **G**      |A|T


### Searching reverse strand through K-mers

**Reference Sequence**: `ATCGGCA`
**K-mer**: `CCG`

Searching for a K-mer on the reverse strand is way simpler. We just reverse complement the k-mer sequence and search for it on the existing strand

**Reverse complement K-mer**: `CGG`

Which is found at genomic coordinates 3,5:  

| 1 | 2 | 3 | 4   | 5 |6|7|
|----------|----------|---------------|----------|--------|-|-|
| A        | T        | **C**             | **G**    | **G**      |C|A

### Searching reverse strand through PWMs (JASPAR format)

**PWM**: 

| Position   | 1   | 2   | 3   | 4   |
|------------|-----|-----|-----|-----|
| A          | 0.2 | 0.1 | 0.4 | 0.3 |
| C          | 0.3 | 0.5 | 0.1 | 0.2 |
| G          | 0.4 | 0.3 | 0.4 | 0.5 |
| T          | 0.1 | 0.1 | 0.1 | 0.0 |

To search for a motif on the reverse strand, the reverse complement of the PWM must be calculated. The reverse complement of the above PWM will be represented by reversing the order of the positions and replacing each nucleotide with its complement. Therefore, the reverse PWM becomes:


Reverse PWM:
| Position   | 4   | 3   | 2   | 1   |
|------------|-----|-----|-----|-----|
| A          | 0.0 | 0.1 | 0.5 | 0.2 |
| C          | 0.2 | 0.1 | 0.5 | 0.3 |
| G          | 0.5 | 0.4 | 0.3 | 0.4 |
| T          | 0.3 | 0.4 | 0.1 | 0.1 |

And the standard process of PWM search is followed afterwards. 

### Summary

MAFin computes conservation by comparing aligned sequences while respecting gaps to maintain the MAF file alignment. The similarity vector, which has the same length as the motif, represents matches and mismatches, and the conservation percentage provides an overall measure of how conserved the motif is across genomes. MAFin also outputs genomic coordinates for each motif in the reference and aligned genomes, allowing for easy tracking of motifs across genomes.


## License
This project is licensed under the [GNU GPL v3](LICENSE).

### Contact

For any questions or support, please contact 
- izg5139@psu.edu
- mpp5977@psu.edu 
- kap6605@psu.edu
- ioannis.mouratidis@psu.edu

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
