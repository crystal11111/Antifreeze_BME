import subprocess
import os
from Bio import SeqIO, AlignIO

afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']

for afp_type in afp_types:
    input_file = f"data/positives/pblast_filtered/pblast_{afp_type}_filtered.fasta"
    
    if not os.path.exists(input_file):
        print(f"Skipping {afp_type}: file not found.")
        continue
        
    # Check if we have enough sequences
    seq_count = len(list(SeqIO.parse(input_file, "fasta")))
    if seq_count < 2:
        print(f"Skipping {afp_type}: only {seq_count} sequence(s)")
        continue
        
    print(f"Processing {afp_type} ({seq_count} sequences)...")
    
    # CD-HIT clustering
    clustered_file = f"data/positives/pblast_filtered/pblast_{afp_type}_c90.faa"
    subprocess.run([
        "cd-hit", "-i", input_file, "-o", clustered_file, 
        "-c", "0.90", "-n", "5"
    ], check=True)
    
    # MSA with ClustalW
    subprocess.run([
        "clustalw", f"-INFILE={clustered_file}"
    ], check=True)
    
    # Check alignment
    try:
        alignment = AlignIO.read(f"{clustered_file.replace('.faa', '.aln')}", "clustal")
        print(f"{afp_type}: {len(alignment)} sequences, length {alignment.get_alignment_length()}")
    except:
        print(f"{afp_type}: alignment failed")

print("Processing complete!")