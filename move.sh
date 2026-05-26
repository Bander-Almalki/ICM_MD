# 1. Target subfolder name
Folder="Data/train/"
SUBFOLDER="Data/train/sequences_fasta/"

# 2. Loop through every file in the current directory
for file in $Folder*; do
    
    # Check if it's an actual file (not a folder) and if it ends with .fasta
    if [ -f "$file" ] && [[ "$file" == *.fasta ]]; then
        echo "Moving: $file -> $SUBFOLDER/"
        mv "$file" "$SUBFOLDER/"
    fi
    
done
