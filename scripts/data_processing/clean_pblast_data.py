import subprocess
import os
from Bio import SeqIO

def normalize_id(seq_id):
    # Strip version suffixes like .1, .2 from FASTA IDs to ensure matching
    # Also stip the id from any pipe characters if present
    if '|' in seq_id:
        seq_id = seq_id.split('|')[1]
    return seq_id.split('.')[0]

afp_types = ['type1', 'type2', 'type3', 'type4']  # AFGP removed

# Load reviewed IDs so they can be excluded from non-reviewed set
reviewed_data_file = "data/raw/uniprotkb_fishAFP.fasta"
with open(reviewed_data_file, "r") as rev_file:
    reviewed_seqs = {normalize_id(record.id) for record in SeqIO.parse(rev_file, "fasta")}

# -----------------------------------------------------------------
# Process each AFP type input file
# -----------------------------------------------------------------
for afp_type in afp_types:
    input_file = f"data/raw/pblast_non_reviewed/pblast_{afp_type}_nonredundant_(nr).fasta"

    if not os.path.exists(input_file):
        print(f"Skipping {afp_type}: file not found.")
        continue

    output_file = f"data/positives/pblast_filtered/pblast_{afp_type}_filtered.fasta"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    filtered_records = []

    for record in SeqIO.parse(input_file, "fasta"):
        desc = record.description.lower()

        # Remove reviewed sequences
        if normalize_id(record.id) in reviewed_seqs:
            print(f"Excluding reviewed sequence: {record.id}")
            continue

        # Remove hypothetical and partial sequences
        if "hypothetical" in desc or "partial" in desc:
            continue

        filtered_records.append(record)

    SeqIO.write(filtered_records, output_file, "fasta")
    print(f"{afp_type}: kept {len(filtered_records)} sequences.")

print("Processing complete!")
