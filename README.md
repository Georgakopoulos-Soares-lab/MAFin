# MAF Motif Search CLI Tool

## Description

This CLI tool allows you to search MAF (Multiple Alignment Format) files for specific motifs using regular expression patterns, K-mers, or Position Weight Matrices (PWMs). It supports searching within reference genomes or across all genomes and can handle both strands if specified. The tool is optimized for performance with support for parallel processing.

## Features

- **Motif Search Options**: Search using regex patterns, K-mers, or PWMs (JASPAR format).
- **Flexible Genome Selection**: Choose to search within the reference genome or across all genomes.
- **Strand Orientation**: Optionally search both strands by considering reverse complements.
- **Parallel Processing**: Utilize multiple processes to speed up the search.
- **Output Formats**: Generate detailed JSON reports of the motifs found.

## Installation

### Prerequisites

- Python 3.6 or higher
- `pip` package manager
-  Biopython ( installed via requirements.txt )

### Using Virtual Environment

1. **Clone the Repository**

   ```bash
   git clone https://github.com/yourusername/maf-motif-search.git
   cd maf-motif-search
    ```

2. **Install inside venv**
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```


## Usage
   ```bash
   python maf_motif_search.py --maf_file path/to/file.maf \
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
                           --processes 8
```

### Searching with K-mers
```bash
python maf_motif_search.py --maf_file data/sample.maf \
                           --kmers data/kmers.txt \
                           --search_in reference \
                           --pvalue_threshold 1e-5 \
                           --processes 2
```

### Searching with PWM (JASPAR format)
```bash
python maf_motif_search.py --maf_file data/sample.maf \
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

Yet to be implemented


### License

We have chosen to apply the Creative Commons 4.0 BY-SA license to Mafin . The description of the license is given below. More information can be found in the Creative Commons web page.


### Contact

For any questions or support, please contact 
- izg5139@psu.edu
- mpp5977@psu.edu 
- kap6605@psu.edu
- ioannis.mouratidis@psu.edu