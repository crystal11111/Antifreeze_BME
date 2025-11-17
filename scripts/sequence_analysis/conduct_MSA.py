from Bio import SeqIO
import subprocess
from Bio import AlignIO

input_faa = "./data/positives/afp_all_c90" # Input file

clustalw_exe = "clustalw2" # Requires being on the environment path

# Parse sequences in file
sequences = list(SeqIO.parse(input_faa+".faa", "fasta"))
print(f"Loaded {len(sequences)} sequences")

# Run ClustalW
subprocess.run(["clustalw2", "-INFILE=" + input_faa + ".faa"], check=True)

# Read the alignment
alignment = AlignIO.read(input_faa+".aln", "clustal")
print(f"Alignment length: {alignment.get_alignment_length()}")
print(alignment)

# Step 5: Save the alignment in another format (optional)
AlignIO.write(alignment, input_faa+"_aligned.faa", "fasta")
