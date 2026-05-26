import mmseqs
import pandas as pd
import argparse

def apply_mmseqs2(sequences, similarity_threshold=0.3):
  """Applies MMseqs2 with a specified similarity threshold to a list of sequences.

  Args:
    sequences: A list of protein sequences.
    similarity_threshold: The desired similarity threshold (0-1).

  Returns:
    A list of clustered sequences.
  """

  # Create a MMseqs2 environment
  env = mmseqs.Environment()

  # Create a database from the sequences
  db = mmseqs.create(sequences, env=env)

  # Perform clustering
  clusters = mmseqs.cluster(sequences, db, env=env, similarity=similarity_threshold, min_seqid=similarity_threshold)

  return clusters

# Example usage
#sequences = ["seq"]  # Your list of sequences
#df = pd.DataFrame({"sequences": sequences})
parser = argparse.ArgumentParser(description="cluster sequences using MMseqs2")
parser.add_argument('--input', required=True, help="Path to the input PDB file")
args = parser.parse_args()
df = pd.read_csv(args.input)

#df=pd.read_csv("/home/light/Downloads/mmseqs2/seq_clusters/TMDOCK_Dimer_Sequences.csv")

clusters = apply_mmseqs2(df["seq"].tolist(), similarity_threshold=0.3)

print(clusters)
