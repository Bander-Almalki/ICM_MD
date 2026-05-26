#!/bin/bash
# Exit immediately if any script fails
set -e

echo "Starting pipeline..."

pdb_dir=Data/train/FMAP_BITOPIC30
fasta_file="Data/train/TMDOCK_Dimer_Sequences.fasta"
database="/home/light/remote_uniref/UniRef30_2023_02"
hhblit_output_dir="Features/train/msa_outputs/"
mmseqs2_output_dir="Features/train/mmseqs2_outputs/"

CCMpred_path="CCMpred/CCMpred/bin/ccmpred"
MSAs_aln="Features/train/dimer_MSAs_aln/"
# new_MSAs_aln="NEW_dimer_MSAs_aln/"
CCMpred_output_dir="CCMpred_output/"
MSAT_output_dir="Features/train/MSAT_output/"
motifs_output_dir="Features/train/motifs_outputs/"
data_fasta="Data/train/sequences_fasta"
pssm_output_dir="Features/train/pssm_outputs/"
lips_output_dir="Features/train/lips_outputs/"
monomer_features_output_dir="Features/train/monomer_features_outputs/"

# 1. Run scripts as standalone command line utilities
# mmseqs2 clustering
#python Scrips/try_clus.py --input $fasta_file $mmseqs2_output_dir
#python Scrips/hhblit_search.py --input $fasta_file $database $hhblit_output_dir
#python Scrips/a3m_to_aln.py --input $hhblit_output_dir $MSAs_aln
#python Scrips/CCMpred.py --input $CCMpred_path $MSAs_aln $CCMpred_output_dir
#python Scrips/msatransformer.py --input $hhblit_output_dir $MSAT_output_dir
#python Scrips/gxxxg_motif.py --input $data_fasta $motifs_output_dir
#python Scrips/sxxxs_motif.py --input $data_fasta $motifs_output_dir
#python Scrips/pxxxp_motif.py --input $data_fasta $motifs_output_dir
#python Scrips/lips_calculation.py --input $MSAs_aln $MSAs_aln
#python Scrips/pssm.py --input $hhblit_output_dir $pssm_output_dir
python extract_monomer_intra_dist_map.py --input $pdb_dir $monomer_features_output_dir




# python data_preprocessing.py --input "U1.pdb" --output "clean.csv"
# python distance_mapping.py --input "clean.csv" --output "matrix.npy"
# python model_training.py --matrix "matrix.npy"

echo "Pipeline complete!"
