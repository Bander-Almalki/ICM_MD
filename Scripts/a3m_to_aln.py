from pathlib import Path
import os

# Set the directory containing your a3m files
DATA_DIR = Path("dimer_MSAs/")
output_dir="dimer_MSAs_aln/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# Loop over all files ending in .a3m in the specified directory
for a3m_path in DATA_DIR.glob("*.a2m"):
    # Generate the output path by replacing the extension with .aln
    aln_path = a3m_path.with_suffix(".aln")
    
    print(f"Processing: {a3m_path.name} -> {aln_path.name}")
    
    cleaned_lines = []
    
    # Read the input file line by line
    with open(a3m_path, "r") as infile:
        for line in infile:
            # 1. egrep -v "^>" (Skip lines starting with '>')
            if line.startswith(">"):
                continue
                
            # 2. sed 's/[a-z]//g' (Remove all lowercase letters)
            # We preserve the trailing newline character
            cleaned_line = "".join([char for char in line if not char.islower()])
            
            cleaned_lines.append(cleaned_line)
            
    # Write the processed lines to the new .aln file
    with open(os.path.join(output_dir,aln_path.name), "w") as outfile:
        outfile.writelines(cleaned_lines)

print("All files processed successfully!")
