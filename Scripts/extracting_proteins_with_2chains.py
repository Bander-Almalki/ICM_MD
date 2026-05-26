'''
This python code takes the TMDOCK predicted homodimer structures and output proteins with 2chains
'''




from Bio.PDB import PDBParser
import os
import pandas as pd

def count_chains_biopython(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    chains = set()
    for model in structure:
        for chain in model:
            chains.add(chain.id)
    if " " in chains:
        chains.remove(" ")
    return len(chains), sorted(chains)



folder="TMDOCK_dimers_pdbs/TMDOCK_DIMER"
c=0
list_of_pdbs_with_two_chains=[]
for filename in os.listdir(folder):
    if filename.endswith(".pdb"):
        number_of_chains, chain_ids = count_chains_biopython(os.path.join(folder, filename))
        if number_of_chains == 2:
            print(f"{filename} has {number_of_chains} chains")
            c+=1
            list_of_pdbs_with_two_chains.append(filename)
print("total", c)
print(len(list_of_pdbs_with_two_chains))

df=pd.DataFrame(list_of_pdbs_with_two_chains,columns=['Protein_Name'])
df.to_csv('Homodimers2chains.csv',index=False)