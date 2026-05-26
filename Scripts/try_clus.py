import subprocess
import os
import pandas as pd
import argparse

def run_mmseqs_clustering(input_fasta, output_prefix, threshold=0.3, coverage=0.8):
    tmp_dir = "tmp_mmseqs"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Construct the MMseqs2 command
    cmd = [
        "mmseqs", "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(threshold),
       
        
    ]

    try:
        print(f"Running MMseqs2 clustering at {int(threshold*100)}% identity...")
        subprocess.run(cmd, check=True)
        print("Clustering complete.")
        
        # Load the TSV output into a pandas DataFrame
        tsv_file = f"{output_prefix}_cluster.tsv"
        df = pd.read_csv(tsv_file, sep='\t', names=['representative', 'member'])
        return df

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        return None

# Usage

parser = argparse.ArgumentParser(description="cluster sequences using MMseqs2")
parser.add_argument('--input',type=str, required=True, help="Path to the fasta file")
args = parser.parse_args()
sequences = args.input

results_df = run_mmseqs_clustering(sequences, "my_clusters")

if results_df is not None:
    print(results_df.head())
    print(f"Total sequences: {len(results_df)}")
    print(f"Total clusters: {results_df['representative'].nunique()}")
