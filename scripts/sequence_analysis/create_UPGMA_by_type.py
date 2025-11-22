from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import os

# Process each AFP type
afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']

for afp_type in afp_types:
    aln_file = f"./data/positives/by_type/{afp_type}_c90.aln"
    
    if not os.path.exists(aln_file):
        print(f"Skipping {afp_type}: alignment file not found")
        continue
    
    try:
        # Read alignment
        alignment = AlignIO.read(aln_file, "clustal")
        
        if len(alignment) < 2:
            print(f"Skipping {afp_type}: need at least 2 sequences")
            continue
            
        print(f"Processing {afp_type}: {len(alignment)} sequences")
        
        # Compute distances and build tree
        calculator = DistanceCalculator("blosum62")
        distance_matrix = calculator.get_distance(alignment)
        
        constructor = DistanceTreeConstructor()
        upgma_tree = constructor.upgma(distance_matrix)
        
        # Save tree
        nwk_file = f"./data/positives/by_type/{afp_type}_c90.nwk"
        Phylo.write(upgma_tree, nwk_file, "newick")
        print(f"Saved {afp_type} tree to {nwk_file}")
        
    except Exception as e:
        print(f"Error processing {afp_type}: {e}")

print("UPGMA tree generation complete!")