from Bio import AlignIO
import pandas as pd
import numpy as np
import os
from collections import Counter

def calculate_conservation_score(column):
    """Calculate conservation score for an alignment column"""
    # Remove gaps
    residues = [res for res in column if res != '-']
    if not residues:
        return 0
    
    # Count residue frequencies
    counts = Counter(residues)
    total = len(residues)
    
    # Shannon entropy-based conservation
    entropy = 0
    for count in counts.values():
        freq = count / total
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    # Convert to conservation score (0 = variable, 1 = conserved)
    max_entropy = np.log2(min(20, len(set(residues))))  # Max possible entropy
    if max_entropy == 0:
        return 1  # Perfectly conserved
    
    conservation = 1 - (entropy / max_entropy)
    return conservation

def extract_msa_features_for_sequence(seq_record, alignment, seq_index):
    """Extract MSA-based features for a single sequence"""
    features = {}
    
    sequence = str(seq_record.seq)
    aln_length = alignment.get_alignment_length()
    
    # Basic sequence features
    features['seq_length'] = len(sequence.replace('-', ''))  # Length without gaps
    features['gap_count'] = sequence.count('-')
    features['gap_percentage'] = (sequence.count('-') / len(sequence)) * 100
    
    # Conservation features
    conservation_scores = []
    residue_conservations = []
    
    for i in range(aln_length):
        column = alignment[:, i]
        conservation = calculate_conservation_score(column)
        conservation_scores.append(conservation)
        
        # Conservation of this sequence's residue at this position
        if sequence[i] != '-':
            residue_conservations.append(conservation)
    
    # Overall conservation statistics
    features['mean_conservation'] = np.mean(conservation_scores)
    features['std_conservation'] = np.std(conservation_scores)
    features['min_conservation'] = np.min(conservation_scores)
    features['max_conservation'] = np.max(conservation_scores)
    
    # Sequence-specific conservation
    if residue_conservations:
        features['seq_mean_conservation'] = np.mean(residue_conservations)
        features['seq_std_conservation'] = np.std(residue_conservations)
        features['highly_conserved_positions'] = sum(1 for c in residue_conservations if c > 0.8)
        features['variable_positions'] = sum(1 for c in residue_conservations if c < 0.2)
    else:
        features['seq_mean_conservation'] = 0
        features['seq_std_conservation'] = 0
        features['highly_conserved_positions'] = 0
        features['variable_positions'] = 0
    
    # Sequence identity to consensus
    consensus_matches = 0
    for i in range(aln_length):
        if sequence[i] != '-':
            column = alignment[:, i]
            residues = [res for res in column if res != '-']
            if residues:
                most_common = Counter(residues).most_common(1)[0][0]
                if sequence[i] == most_common:
                    consensus_matches += 1
    
    seq_len_no_gaps = len(sequence.replace('-', ''))
    features['consensus_identity'] = (consensus_matches / seq_len_no_gaps * 100) if seq_len_no_gaps > 0 else 0
    
    return features

def main():
    print("=== Extracting MSA Conservation Features ===")
    
    afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']
    all_features = []
    
    for afp_type in afp_types:
        aln_file = f"data/positives/by_type/{afp_type}_c90.aln"
        
        if not os.path.exists(aln_file):
            print(f"Skipping {afp_type}: alignment file not found")
            continue
            
        try:
            print(f"Processing {afp_type}...")
            
            # Read alignment
            alignment = AlignIO.read(aln_file, "clustal")
            
            if len(alignment) < 2:
                print(f"  Skipping {afp_type}: need at least 2 sequences")
                continue
            
            # Extract features for each sequence
            for i, seq_record in enumerate(alignment):
                features = extract_msa_features_for_sequence(seq_record, alignment, i)
                
                # Extract sequence ID from record name
                seq_id = seq_record.id.split('|')[1] if '|' in seq_record.id else seq_record.id
                features['sequence_id'] = seq_id
                features['afp_type'] = afp_type.upper()
                
                all_features.append(features)
                
            print(f"  Extracted features for {len(alignment)} sequences")
            
        except Exception as e:
            print(f"Error processing {afp_type}: {e}")
    
    # Save to CSV
    if all_features:
        df = pd.DataFrame(all_features)
        output_file = "data/features/msa_features.csv"
        os.makedirs("data/features", exist_ok=True)
        df.to_csv(output_file, index=False)
        print(f"\nSaved MSA features to {output_file}")
        print(f"Total sequences: {len(df)}")
        print(f"Features per sequence: {len(df.columns)-2}")  # Exclude seq_id and afp_type
        
        # Show sample
        print("\nSample features:")
        print(df[['sequence_id', 'afp_type', 'seq_length', 'mean_conservation', 'consensus_identity']].head())
    else:
        print("No features extracted!")

if __name__ == "__main__":
    main()