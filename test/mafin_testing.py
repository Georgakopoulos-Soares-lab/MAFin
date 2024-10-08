import os
import subprocess
import json

# Define directory paths
BASE_DIR = "/storage/group/izg5139/default/multiple_alignment/maffin_git"
TEST_DIR = os.path.join(BASE_DIR, "test")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
REF_GENOME_DIR = os.path.join(RESULTS_DIR, "ref_genome")

# Paths for test data
MAF_FILE = os.path.join(TEST_DIR, "test_alignment.maf")
JASPAR_FILE = os.path.join(TEST_DIR, "test_CCATATATAG_pwm.jaspar")
KMERS_FILE = os.path.join(TEST_DIR, "test_16mer.txt")

# Define the path to the mafin.py script
MAFFIN_SCRIPT = os.path.join(BASE_DIR, "mafin.py")

def run_maffin(args):
    """Runs the mafin.py script with the given arguments."""
    command = ["python", MAFFIN_SCRIPT] + args
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"Error running mafin.py: {result.stderr.decode()}")
        raise RuntimeError("mafin.py failed")
    return result.stdout.decode()

def test_pwm_trivial():
    # Run PWM test
    args = ["--maf_file", MAF_FILE, "--jaspar", JASPAR_FILE, "--purge_results_dir"]
    run_maffin(args)

    # Assert the resulting JSON contains 34 entries
    result_json_file = os.path.join(REF_GENOME_DIR, "test_CCATATATAG_pwm_MA0001.3_motif_hits.json")
    with open(result_json_file, 'r') as f:
        data = json.load(f)
        assert len(data) == 34, f"Expected 34 entries, found {len(data)}"


def test_pwm_background_freq_modified():
    # Run PWM test
    args = ["--maf_file", MAF_FILE, "--jaspar", JASPAR_FILE, "--purge_results_dir" , "--background_frequencies", '0.25', '0.25', '0.15', '0.35']
    run_maffin(args)

    # Assert the resulting JSON contains 43 entries
    result_json_file = os.path.join(REF_GENOME_DIR, "test_CCATATATAG_pwm_MA0001.3_motif_hits.json")
    with open(result_json_file, 'r') as f:
        data = json.load(f)
        assert len(data) == 43, f"Expected 43 entries, found {len(data)}"



def test_regex_trivial():
    # Run regex test
    args = ["--maf_file", MAF_FILE, "--regex", "A{5}T{4}(GGG)?C{2}", "--purge_results_dir"]
    run_maffin(args)

    # Assert the resulting JSON contains 6 entries
    result_json_file = os.path.join(REF_GENOME_DIR, "regex_search_A_5_T_4_GGG_C_2__motif_hits.json")
    with open(result_json_file, 'r') as f:
        data = json.load(f)
        assert len(data) == 6, f"Expected 6 entries, found {len(data)}"

def test_kmer_trivial():
    # Run kmer test
    args = ["--maf_file", MAF_FILE, "--kmers", KMERS_FILE, "--purge_results_dir"]
    run_maffin(args)

    # Assert the resulting JSON contains 6 entries
    result_json_file = os.path.join(REF_GENOME_DIR, "test_16mer_test_16mer_motif_hits.json")
    with open(result_json_file, 'r') as f:
        data = json.load(f)
        assert len(data) == 6, f"Expected 6 entries, found {len(data)}"

# Run the tests
if __name__ == "__main__":
    test_pwm_trivial()
    test_regex_trivial()
    test_kmer_trivial()
    test_pwm_background_freq_modified()
    print("All tests passed successfully.")
