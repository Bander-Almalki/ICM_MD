import re
import pandas as pd
import numpy as np
import argparse
import os
from Bio import SeqIO


motif_small = "[GASC].{3}[GASC]"
def find_sxxxs_motif(seq,motif_len):
  match_start_list = []
  match_end_list = [0] * motif_len
  # counter for the matches
  match_number = 0
  result_dict = {}

  for start in range(len(seq) - motif_len):
      # for a SmallxxxSmall motif, the end is 4 residues later
      end = start + motif_len
      # get the matched segment
      segment = seq[start:end + 1]
      # check if the segment contains a motif
      match = re.match(motif_small, segment)
      if match:
          # classify position as start of a motif
          match_start_list.append(1)
          match_end_list.append(1)
      else:
          match_start_list.append(0)
          match_end_list.append(0)
  # add the final zeros on the end of the list, so it's length matches the original sequence
  match_start_list = match_start_list + [0] * motif_len

  match_start_arr = np.array(match_start_list)
  match_end_arr = np.array(match_end_list)

  list_residues_in_motif = match_start_arr + match_end_arr
  list_residues_in_motif[list_residues_in_motif > 1] = 1
  return list_residues_in_motif



parser = argparse.ArgumentParser()
parser.add_argument('--input',nargs=2, required=True)
args = parser.parse_args()

folder_path=args.input[0]
output_dir=args.input[1]
sub_dir="SxxxS/"


if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(output_dir+sub_dir):
    os.makedirs(output_dir+sub_dir)

for filename in os.listdir(folder_path):
    if filename.endswith(".fasta"):
        msa_path = os.path.join(folder_path, filename)
        seq=str(SeqIO.read(msa_path, "fasta").seq)
        output_file=output_dir+sub_dir+filename.split(".")[0]+".csv"
        sxxxs_score=find_sxxxs_motif(seq,5)
        pd.DataFrame(sxxxs_score).to_csv(output_file)
        print(f"{filename} processed")

# seq="GKVIGIAVMALLLASALTLLIGTVSNMFQ"
# find_sxxxs_motif(seq,4)
