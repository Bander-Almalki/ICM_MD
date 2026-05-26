import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB
import os
import argparse


def get_distance_matrix(pdb_file_path):
    # Initialize the parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file_path)
    
    # Extract Ca atoms from the first model, first chain
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Ensure it's a standard amino acid residue (ignores heteroatoms/water)
                if PDB.is_aa(residue):
                    if 'CA' in residue:
                        ca_atoms.append(residue['CA'])
            # Since it's a monomer, we can break after the first chain
            break
        break

    num_residues = len(ca_atoms)
    print(f"Extracted {num_residues} C-alpha atoms.")

    # Initialize an empty matrix
    distance_matrix = np.zeros((num_residues, num_residues))

    # Calculate Euclidean distances
    for i in range(num_residues):
        for j in range(i, num_residues):
            # Biopython overrides the subtraction operator to calculate distance
            distance = ca_atoms[i] - ca_atoms[j]
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  # The matrix is symmetric

    return distance_matrix

def plot_distance_map(matrix):
    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='viridis', origin='lower')
    plt.colorbar(label='Distance (Å)')
    plt.title('Intra-Chain $C_\\alpha$ Distance Map')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()



parser = argparse.ArgumentParser()
parser.add_argument('--input',nargs=2, required=True)
args = parser.parse_args()

folder_path=args.input[0]
output_dir=args.input[1]


if not os.path.exists(output_dir):
    os.makedirs(output_dir)


for filename in os.listdir(folder_path):
    if filename.endswith(".pdb"):
        pdb_file = os.path.join(folder_path, filename)
        output_file = os.path.join(output_dir, filename.split('.')[0]+'_intra_dist_map.csv')
        dist_matrix = get_distance_matrix(pdb_file)
        #plot_distance_map(dist_matrix)
        np.savetxt(output_file, dist_matrix, delimiter=',')
        print(f"Extracted data for {filename} and saved to {output_file}")



# # --- Execution ---
# output_directory="monomer_structure_features/"
# if not os.path.exists(output_directory):
#   os.makedirs(output_directory)

# directory='FMAP_BITOPIC30/'
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         pdb_file = os.path.join(directory, filename)
#         output_file = os.path.join(output_directory, filename.split('.')[0]+'_intra_dist_map.csv')
#         dist_matrix = get_distance_matrix(pdb_file)
#         #plot_distance_map(dist_matrix)
#         np.savetxt(output_file, dist_matrix, delimiter=',')
#         print(f"Extracted data for {filename} and saved to {output_file}")
#     break



