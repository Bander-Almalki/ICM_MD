import torch
import esm
from Bio import SeqIO
import itertools
import os 
import pickle
import argparse

def clean_msa_sequence(sequence):
    """
    Removes lowercase letters (insertions in A3M) and replaces 
    dots with dashes to ensure a clean alignment for the transformer.
    """
    # 1. Remove lowercase letters (A3M insertion format)
    cleaned = "".join([c for c in sequence if not c.islower()])
    # 2. Standardize gaps
    return cleaned.replace(".", "-")

def process_msa_for_transformer(msa_path, num_sequences=128):
    msa_data = []
    
    # Load sequences
    records = list(itertools.islice(SeqIO.parse(msa_path, "fasta"), num_sequences))
    
    # Determine target length from the first sequence (the query)
    query_seq = clean_msa_sequence(str(records[0].seq))
    target_len = len(query_seq)
    
    for record in records:
        cleaned = clean_msa_sequence(str(record.seq))
        
        # Validation: If a sequence still doesn't match the query length, 
        # it's usually because the A3M was malformed or it's a partial hit.
        if len(cleaned) != target_len:
            # Option: Force padding/truncating if strictly necessary, 
            # though usually cleaning lowercase is enough.
            cleaned = cleaned[:target_len].ljust(target_len, "-")
            
        msa_data.append((record.id, cleaned))
    
    return msa_data

def get_msa_row_attention(msa_path, num_sequences=64):
    # 1. Load model and alphabet
    model, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    model = model.eval()
    batch_converter = alphabet.get_batch_converter()

    if torch.cuda.is_available():
        model = model.cuda()

    # 2. Prepare MSA data (subsampling is critical for memory)

    # msa_data = []
    # for record in itertools.islice(SeqIO.parse(msa_path, "fasta"), num_sequences):
    #     # ESM is picky about gap characters; ensure they are '-'
    #     sequence = str(record.seq).replace(".", "-")
    #     msa_data.append((record.id, sequence))
    msa_data = process_msa_for_transformer(msa_path)
    #query_sequence = msa_data[0] # Grab just the first (id, seq) tuple

    # 2. IMPORTANT: The model expects a list containing the MSA list
    # This creates a 'batch' where the batch size is 1 MSA.
    batch_data = [msa_data]

    # 3. Call the converter
    # If it still fails, try removing the labels from the unpack:
    # results = batch_converter(batch_data)
    # batch_labels, batch_strs, batch_tokens = results
    batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)

    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()

    # 3. Forward pass with attention weights enabled
    with torch.no_grad():
        results = model(batch_tokens, return_contacts=True, need_head_weights=True)
    
    print(results.keys())
    # 4. Extract Attention
    # Shape: [layers, batch, heads, num_sequences, seq_len, seq_len]
    attentions = results["row_attentions"]
    last_layer_attn = attentions[-1,:,:,1:,1:]
    print(f"Attention Shape: {last_layer_attn.shape}")
    
    # The MSA Transformer alternates between Row and Column attention.
    # Typically, even layers (0, 2, 4...) contain row attention.
    # Let's extract row attention from the last layer (layer 11)
    row_attention = last_layer_attn

    return row_attention

row_att={}

parser = argparse.ArgumentParser(description="Extract MSA Transformer features")
parser.add_argument('--input',nargs=2, required=True)
args = parser.parse_args()

folder=args.input[0]
output_dir=args.input[1]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#itrate over msa files in a folder
#folder ="../dimer_MSAs/"
for filename in os.listdir(folder):
    if filename.endswith(".a2m"):
        msa_path = os.path.join(folder, filename)
        row_attn=get_msa_row_attention(msa_path)
        L=row_attn.shape[2]
        reshaped_attn = row_attn.permute(2, 3, 0, 1)
        final_tensor = reshaped_attn.reshape(L, L, 144)

        print(f"filename: {filename}")
        print(f"Row Attention Shape: {final_tensor.shape}")
        #store the tensor into a dictionary 
        row_att[filename]=final_tensor

#save the dictionary as a pickle file
with open(os.path.join(output_dir, "msat_row_att.pickle"), "wb") as f: 
    pickle.dump(row_att, f)


