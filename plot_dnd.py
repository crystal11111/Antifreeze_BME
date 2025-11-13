from Bio import Phylo
import matplotlib.pyplot as plt

# Path to .dnd
dnd_file = "./data/positives/afp_all_c90.dnd" 

# ClustalW produces Newick format of tree
tree = Phylo.read(dnd_file, "newick")

# Plot the tree graphically
plt.figure(figsize=(8, 6))
Phylo.draw(tree, do_show=False)

plt.title("ClustalW Guide Tree (.dnd)", fontsize=14)
plt.tight_layout()
plt.show()
