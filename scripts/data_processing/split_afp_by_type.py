import pandas as pd
from Bio import SeqIO
import os

# Read metadata
meta_df = pd.read_csv("data/meta/uniprotkb_fishAFP.tsv", sep="\t")

# Create type mapping from protein names
type_mapping = {}
for _, row in meta_df.iterrows():
    entry_id = row['Entry']
    protein_name = row['Protein names']
    
    # Determine AFP type from protein name
    if 'Type-2' in protein_name or 'Type II' in protein_name:
        afp_type = 'Type2'
    elif 'Type-3' in protein_name or 'Type III' in protein_name:
        afp_type = 'Type3'
    elif 'Type-4' in protein_name or 'Type IV' in protein_name:
        afp_type = 'Type4'
    elif 'glycoprotein' in protein_name or 'AFGP' in protein_name:
        afp_type = 'AFGP'
    else:
        afp_type = 'Type1'  # Default for others
    
    type_mapping[entry_id] = afp_type

# Create output directory
os.makedirs("data/positives/by_type", exist_ok=True)

# Initialize file handles
type_files = {}
for afp_type in set(type_mapping.values()):
    type_files[afp_type] = open(f"data/positives/by_type/{afp_type.lower()}_fish.fasta", "w")

# Parse FASTA and split by type
for record in SeqIO.parse("data/raw/uniprotkb_fishAFP.fasta", "fasta"):
    entry_id = record.id.split("|")[1]  # Extract UniProt ID
    afp_type = type_mapping.get(entry_id, 'Type1')
    SeqIO.write(record, type_files[afp_type], "fasta")

# Close files
for f in type_files.values():
    f.close()

print("AFP sequences split by type:")
for afp_type in type_files.keys():
    count = len(list(SeqIO.parse(f"data/positives/by_type/{afp_type.lower()}_fish.fasta", "fasta")))
    print(f"{afp_type}: {count} sequences")