from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

file = "./data/positives/afp_all_c90" # Expects a file with an ending of .aln to exist

# Read alignment generated from clustalw
alignment = AlignIO.read(file+".aln", "clustal")

# Compute pairwise distances
calculator = DistanceCalculator("blosum62")
distance_matrix = calculator.get_distance(alignment)
print(distance_matrix)

# Build UPGMA tree
constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(distance_matrix)

# Plot tree
Phylo.draw(upgma_tree, do_show=False)
plt.title("UPGMA Tree from ClustalW Alignment", fontsize=14)
plt.tight_layout()
plt.show()

# Save Tree
Phylo.write(upgma_tree, file + ".nwk", "newick")
