import os
import re
import pandas as pd
import numpy as np
import argparse

def find_gxxxg_motif(seq,motif_len):
  motif_ss = "G...G"
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
      match = re.match(motif_ss, segment)
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


print(find_gxxxg_motif("SSSKSNLGLIVGSAIGSLLAVVFLGSCFVLYKK",5))
