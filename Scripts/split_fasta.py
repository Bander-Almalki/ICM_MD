from pathlib import Path
import re

# Set your input file name and output folder
INPUT_FASTA = "TMDOCK_Dimer_Sequences.fasta"
OUTPUT_DIR = Path("Data/train")

# Create the output folder if it doesn't exist yet
OUTPUT_DIR.mkdir(exist_ok=True)

current_header = None
current_seq_lines = []

def save_current_sequence(header, lines):
    if not header or not lines:
        return
    
    # Sanitize the header to create a safe file name
    # Strips the '>' and replaces spaces/special characters with underscores
    clean_name = re.sub(r'[^\w\-_]', '_', header[1:].split()[0])
    output_path = OUTPUT_DIR / f"{clean_name}.fasta"
    
    with open(output_path, "w") as out_file:
        out_file.write(header + "\n")
        out_file.writelines(lines)

# Read the big file and split it line by line
with open(INPUT_FASTA, "r") as infile:
    for line in infile:
        line_stripped = line.strip()
        if not line_stripped:
            continue  # Skip empty lines
            
        if line_stripped.startswith(">"):
            # If we already have a sequence collected, save it before moving on
            save_current_sequence(current_header, current_seq_lines)
            
            # Reset trackers for the next sequence
            current_header = line_stripped
            current_seq_lines = []
        else:
            # Append sequence lines (keeping their original formatting/newlines)
            current_seq_lines.append(line)

    # Don't forget to save the very last sequence in the file!
    save_current_sequence(current_header, current_seq_lines)

print(f"Done! Separate files have been saved into the '{OUTPUT_DIR}' folder.")
