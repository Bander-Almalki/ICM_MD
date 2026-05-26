import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description="Extract CCMpred features")
parser.add_argument('--input',nargs=3, required=True, help="Path to CCMpred and aln files")
args = parser.parse_args()

CCMpred_path = args.input[0]
aln_dir=args.input[1]
output_dir=args.input[2]

# CCMpred_path = "CCMpred/CCMpred/bin/ccmpred"
# aln_dir="dimer_MSAs_aln/"
# output_dir="CMpred_outputs"
#file=os.path.join(aln_dir,"COX8_YEAST.fasta.aln")


if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file in os.listdir(aln_dir):
    if file.endswith(".aln"):
        aln_path = os.path.join(aln_dir, file)
        output_file = os.path.join(output_dir, file.replace(".aln", ".mat"))
        subprocess.run([CCMpred_path, aln_path, output_file])

