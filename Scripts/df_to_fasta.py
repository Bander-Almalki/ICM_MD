import pandas as pd
import os


df = pd.read_csv("TMDOCK_Dimer_Sequences.csv")
fasta_path = "TMDOCK_Dimer_Sequences2.fasta"
with open(fasta_path, "w") as f:
    for _, row in df.iterrows():
        f.write(f">{row['prot_name']}\n{row['seq']}\n")

print(f"FASTA file saved to {fasta_path}")
