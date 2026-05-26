from Bio import SeqIO
import os
import pandas as pd
import numpy as np
import pickle
import argparse


def generate_pssm(fasta_file):
    # Load MSA from file
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))

    # Get sequence length
    seq_len = len(sequences[0].seq)

    # Initialize PSSM matrix
    pssm = np.zeros((20, seq_len))

    # Amino acid order
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'

    # Iterate over positions
    for i in range(seq_len):
        # Get amino acid frequencies at current position
        freqs = {}
        for seq in sequences:
            amino_acid = seq.seq[i]
            if amino_acid != '-':  # Ignore gaps
                freqs[amino_acid] = freqs.get(amino_acid, 0) + 1

        # Normalize frequencies
        total = sum(freqs.values())
        for amino_acid in amino_acids:
            pssm[amino_acids.index(amino_acid), i] = freqs.get(amino_acid, 0) / total

    return pssm


def print_pssm(pssm):
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'
    for i, row in enumerate(pssm):
        print(f"{amino_acids[i]}: {' '.join(f'{x:.2f}' for x in row)}")



parser = argparse.ArgumentParser()
parser.add_argument('--input',nargs=2, required=True)
args = parser.parse_args()

folder_path=args.input[0]
output_dir=args.input[1]


if not os.path.exists(output_dir):
    os.makedirs(output_dir)


for filename in os.listdir(folder_path):
    if filename.endswith(".a2m"):
        msa_path = os.path.join(folder_path, filename)
        output_file=output_dir+filename.split(".")[0]+".csv"
        pssm = generate_pssm(msa_path)
        pssm = pssm.T
        pd.DataFrame(pssm).to_csv(output_file)
        print(f"{filename} processed")






# output_dir="PSSM_Feature/"
# if not os.path.exists(output_dir):
# 	os.makedirs(output_dir)


# directory='dimer_MSAs/'
# pssm_res={}

# for filename in os.listdir(directory):
# 	if filename.endswith("fasta.a2m"):
# 		fasta_path = os.path.join(directory, filename)
# 		output_file=output_dir+filename.split(".")[0]+".csv"
# 		pssm = generate_pssm(fasta_path)
# 		pssm = pssm.T
# 		pd.DataFrame(pssm).to_csv(output_file)
# 		pssm_res[filename.split(".")[0]]=pssm
# 		print(f"{filename} processed")
	

# with open (f'{output_dir}pssm_res.pickle', 'wb') as handle:
# 	pickle.dump(pssm_res, handle, protocol=pickle.HIGHEST_PROTOCOL)