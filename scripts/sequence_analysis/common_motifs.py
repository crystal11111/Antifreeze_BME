import subprocess
import os
from collections import Counter
from Bio import SeqIO, AlignIO



def extract_motifs(alignment, threshold):
    """
    Extracts stretches of residues that meet the conservation threshold.
    
    Args:
        alignment (Bio.Align.MultipleSeqAlignment): The MSA object.
        threshold (float): The minimum ratio of the most common residue 
                           (excluding gaps) for a column to be conserved.
                           
    Returns:
        list: A list of the extracted motif strings.
    """
    alignment_length = alignment.get_alignment_length()
    motifs = []
    current_motif = ""

    for i in range(alignment_length):
        column = alignment[:, i] # Get the residues at this position
        non_gap_residues = [res for res in column if res != '-']
        if not non_gap_residues:
            # If the entire column is gaps, treat it as non-conserved
            if current_motif:
                motifs.append(current_motif)
            current_motif = ""
            continue
        residue_counts = Counter(non_gap_residues)
        most_common_residue_count = 0
        if residue_counts:
            most_common_residue_count = residue_counts.most_common(1)[0][1]
    
        conservation_ratio = most_common_residue_count / len(non_gap_residues)

        if conservation_ratio >= threshold:
            representative_residue = residue_counts.most_common(1)[0][0]
            current_motif += representative_residue
        else:
            if current_motif:
                motifs.append(current_motif)
            current_motif = ""
    if current_motif:
        motifs.append(current_motif)
        
    # Filter out short motifs (e.g., length 3 or less)
    return [m for m in motifs if len(m) > 3]


if __name__ == "__main__":

    afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']
    conservation_threshold = 0.90 

    for afp_type in afp_types:
        input_file = f"../../data/positives/by_type/{afp_type}_fish.fasta"
        alignment_file = f"../../data/positives/by_type/{afp_type}_c90.aln"

        output_content = []
        if not os.path.exists(input_file):
            print("path doesnt exist")
            continue
        # Check if we have enough sequences
        try:
            seq_count = len(list(SeqIO.parse(input_file, "fasta")))
        except:
            seq_count = 0
            
        if seq_count < 2:
            print(f"Skipping {afp_type}: only {seq_count} sequence(s)")
            continue
            
        print(f"\nProcessing {afp_type} ({seq_count} sequences)...")
        if not os.path.exists(alignment_file):
            print(f"{afp_type}: Alignment file (.aln) not found. Run CD-HIT and ClustalW first.")
            continue
        
        alignment = AlignIO.read(alignment_file, "clustal")
        print(f"-> Alignment: {len(alignment)} sequences, length {alignment.get_alignment_length()}")
        
        motifs = extract_motifs(alignment, conservation_threshold)
        

        print(f"-> Extracted Motifs (Conservation >= {int(conservation_threshold*100)}%):")
        if motifs:
            for motif in motifs:
                print(f"   - {motif} (Length: {len(motif)})")
                output_content.append(motif)
                output_content.append("\n" + "="*50)
        else:
            print("   - No motifs found matching the threshold.")
                
        
        output_folder = rf"..\..\data\motifs"
        output_path = os.path.join(output_folder, f"{afp_type}_motifs.txt")

        with open(output_path, 'w') as f:
            f.write('\n'.join(output_content) + '\n')

        print("\nProcessing complete!")