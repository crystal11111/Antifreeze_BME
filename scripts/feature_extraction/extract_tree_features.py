from Bio import Phylo
import pandas as pd
import numpy as np
import os

def extract_tree_distances(tree_file):
    """Extract pairwise distances from UPGMA tree"""
    tree = Phylo.read(tree_file, "newick")
    terminals = tree.get_terminals()
    
    distances = {}
    for i, term1 in enumerate(terminals):
        for j, term2 in enumerate(terminals):
            if i <= j:  # Only upper triangle + diagonal
                dist = tree.distance(term1, term2)
                seq1_id = term1.name.split('|')[1] if '|' in term1.name else term1.name
                seq2_id = term2.name.split('|')[1] if '|' in term2.name else term2.name
                distances[(seq1_id, seq2_id)] = dist
    
    return distances

def extract_tree_features_for_sequence(seq_id, distances):
    """Extract tree-based features for a single sequence"""
    features = {}
    
    # Get all distances involving this sequence
    seq_distances = []
    for (id1, id2), dist in distances.items():
        if id1 == seq_id or id2 == seq_id:
            if id1 != id2:  # Exclude self-distance (0)
                seq_distances.append(dist)
    
    if seq_distances:
        features[f'tree_min_dist'] = min(seq_distances)
        features[f'tree_max_dist'] = max(seq_distances)
        features[f'tree_mean_dist'] = np.mean(seq_distances)
        features[f'tree_std_dist'] = np.std(seq_distances)
    else:
        # Single sequence in tree
        features[f'tree_min_dist'] = 0
        features[f'tree_max_dist'] = 0
        features[f'tree_mean_dist'] = 0
        features[f'tree_std_dist'] = 0
    
    return features

def main():
    print("=== Extracting UPGMA Tree Features ===")
    
    afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']
    all_features = []
    
    for afp_type in afp_types:
        tree_file = f"data/positives/by_type/{afp_type}_c90.nwk"
        
        if not os.path.exists(tree_file):
            print(f"Skipping {afp_type}: tree file not found")
            continue
            
        try:
            print(f"Processing {afp_type}...")
            
            # Extract distances
            distances = extract_tree_distances(tree_file)
            
            # Get unique sequence IDs
            seq_ids = set()
            for (id1, id2) in distances.keys():
                seq_ids.add(id1)
                seq_ids.add(id2)
            
            # Extract features for each sequence
            for seq_id in seq_ids:
                features = extract_tree_features_for_sequence(seq_id, distances)
                features['sequence_id'] = seq_id
                features['afp_type'] = afp_type.upper()
                all_features.append(features)
                
            print(f"  Extracted features for {len(seq_ids)} sequences")
            
        except Exception as e:
            print(f"Error processing {afp_type}: {e}")
    
    # Save to CSV
    if all_features:
        df = pd.DataFrame(all_features)
        output_file = "data/features/tree_features.csv"
        os.makedirs("data/features", exist_ok=True)
        df.to_csv(output_file, index=False)
        print(f"\nSaved tree features to {output_file}")
        print(f"Total sequences: {len(df)}")
        print(f"Features per sequence: {len(df.columns)-2}")  # Exclude seq_id and afp_type
        
        # Show sample
        print("\nSample features:")
        print(df.head())
    else:
        print("No features extracted!")

if __name__ == "__main__":
    main()