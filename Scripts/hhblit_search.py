import os
import subprocess
from Bio import SeqIO  # Requires: pip install biopython
import argparse


def run_hhblits_batch(input_fasta, db_path, output_dir="msa_outputs"):
    """
    Runs HHblits on each sequence in a FASTA file.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through each sequence in the FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_id = record.id
        # Create a temporary single-sequence FASTA file for HHblits
        temp_query = f"temp_{seq_id}.fasta"
        with open(temp_query, "w") as f:
            f.write(f">{seq_id}\n{record.seq}\n")
        
        # Define output paths
        a3m_output = os.path.join(output_dir, f"{seq_id}.a3m")
        hhr_output = os.path.join(output_dir, f"{seq_id}.hhr")

        # HHblits command with requested parameters
        # -n: iterations, -e: E-value, -id: max sequence identity
        cmd = [
            "hhblits",
            "-i", temp_query,
            "-d", db_path,
            "-oa3m", a3m_output,
            "-o", hhr_output,
            "-n", "3",
            "-e", "0.001",
            "-id", "90",
            "-cpu", "4"  # Adjust based on your available cores
        ]

        try:
            print(f"Processing sequence: {seq_id}...")
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {seq_id}: {e}")
        finally:
            # Clean up temporary query file
            if os.path.exists(temp_query):
                os.remove(temp_query)

    print("\nBatch processing complete.")

# Configuration
parser = argparse.ArgumentParser(description="Extract MultiMonomer sequences alignments using hhblits")
parser.add_argument('--input',nargs=2, required=True, help="Path to the fasta file and database")
args = parser.parse_args()
sequences = args.input[0]
database = args.input[1]
output_dir = args.input[2]

run_hhblits_batch(sequences, database, output_dir)


#INPUT_FILE = "my_clusters_rep_seq.fasta"
#DATABASE_PATH = "/home/light/remote_uniref/UniRef30_2023_02" # Change to your DB path

run_hhblits_batch(INPUT_FILE, DATABASE_PATH)
