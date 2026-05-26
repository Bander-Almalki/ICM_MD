import os 
import argparse

def convert_to_ccmpred_format(input_a3m, output_fasta):
    records = []
    with open(input_a3m, "r") as f:
        lines = f.readlines()

    curr_id, curr_seq = "", ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if curr_id:
                # CCMpred requires: No lowercase, dots to dashes
                clean = "".join([c for c in curr_seq if not c.islower()]).replace(".", "-")
                records.append(f"{clean}\n")
            curr_id = line[1:]
            curr_seq = ""
        else:
            curr_seq += line
    
    if curr_id:
        clean = "".join([c for c in curr_seq if not c.islower()]).replace(".", "-")
        records.append(f">{curr_id}\n{clean}\n")

    with open(output_fasta, "w") as f:
        f.writelines(records)


#iterate through all a3m files in a folder
parser = argparse.ArgumentParser(description="a3m to aln converter")
parser.add_argument('--input',nargs=2, required=True)
args = parser.parse_args()
folder = args.input[0]
output_folder = args.input[1]



#folder = "msa_outputs"
#output_folder = "NEW_dimer_MSAs_aln/"


if not os.path.exists(output_folder):
    os.makedirs(output_folder)
for filename in os.listdir(folder):
    if filename.endswith(".a2m"):
        input_a3m = os.path.join(folder, filename)
        output_fasta = os.path.join(output_folder, filename.replace(".a2m", ".aln"))
        convert_to_ccmpred_format(input_a3m, output_fasta)
        print(f"File {input_a3m} converted to {output_fasta}")


