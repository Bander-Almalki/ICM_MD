import numpy as np
import pandas as pd
import os

def convert_distance_map_to_pairs(distance_matrix, residue_labels=None,fpath=None):
    
    print(distance_matrix.shape)
    output_csv=fpath.split('.')[0]+'res_pair_intra_dist.csv'
    """
    Converts an L x L distance matrix into a pairwise feature list.
    
    Parameters:
    - distance_matrix: np.ndarray of shape (L, L)
    - residue_labels: list of strings of length L (optional, e.g., ['ASN_697', 'PRO_698', ...])
    - output_csv: string, filename for the saved CSV
    """
    L = distance_matrix.shape[0]
    
    # Get the row and column indices for the unique upper-triangle pairs
    # k=1 shifts above the diagonal to exclude self-distances (i == j)
    row_idx, col_idx = np.triu_indices(L, k=1)
    
    # Extract distances for these unique pairs
    pair_distances = distance_matrix[row_idx, col_idx]
    
    # Initialize the feature dictionary
    features = {
        'Residue1_Index': row_idx,
        'Residue2_Index': col_idx
    }
    
    # Add mapped labels if provided (e.g. "ALA_701")
    if residue_labels is not None and len(residue_labels) == L:
        features['Residue1_Label'] = [residue_labels[i] for i in row_idx]
        features['Residue2_Label'] = [residue_labels[j] for j in col_idx]
        
    features['Distance_A'] = pair_distances
    
    # Convert to DataFrame and export to CSV
    df = pd.DataFrame(features)
    df.to_csv(output_csv, index=False)
    
    print(f"Successfully flattened {len(df)} unique pairs into '{output_csv}'.")
    return df

# --- Example Usage ---
# # Assuming 'distance_map' is your already extracted L x L numpy array:
# df_features = convert_distance_map_to_pairs(distance_map, residue_labels=None)

dir="monomer_structure_features/"
output_dir="Res_pair_intra_distances"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file in os.listdir(dir):
    if file.endswith("ntra_dist_map.csv"):
        file_path = os.path.join(dir, file)
        distance_matrix = np.genfromtxt(file_path, delimiter=',')
        df_features = convert_distance_map_to_pairs(distance_matrix, residue_labels=None,fpath=file_path)
        output_file = os.path.join(output_dir, file)
        df_features.to_csv(output_file, index=False)