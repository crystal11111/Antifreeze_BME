from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import os
from medoid_tools import get_medoids_for_types

'''
- Loads per-type UPGMA tree medoids (Type1, Type2, Type3, Type4, AFGP) from positives.
- For every sequence (AFP + NON_AFP), globally aligns the sequence to each type’s medoid
  and emits per-type alignment features.
'''

def align_to_medoid(seq, medoid_seq, aligner):
    """Compute alignment stats to a medoid sequence."""
    seq_u = seq.upper()
    med_u = medoid_seq.upper()
    aln = aligner.align(seq_u, med_u)[0]
    q_blocks, r_blocks = aln.aligned
    
    matches = mismatches = 0
    for (qs, qe), (rs, re) in zip(q_blocks, r_blocks):
        for a, b in zip(seq_u[qs:qe], med_u[rs:re]):
            if a == b:
                matches += 1
            else:
                mismatches += 1
    
    q_aligned = sum(qe - qs for qs, qe in q_blocks)
    r_aligned = sum(re - rs for rs, re in r_blocks)
    aligned_len = min(q_aligned, r_aligned)
    
    if aligned_len == 0:
        return {'score': aln.score, 'identity': 0.0, 'coverage': 0.0, 'gap_rate': 1.0}
    
    ins_q = max(0, q_aligned - aligned_len)
    ins_r = max(0, r_aligned - aligned_len)
    gaps = ins_q + ins_r
    
    return {
        'score': aln.score,
        'identity': matches / aligned_len,
        'coverage': aligned_len / max(1, len(seq_u)),
        'gap_rate': gaps / max(1, aligned_len + gaps)
    }

def extract_medoid_features(sequence, seq_id, medoid_info, label):
    """Extract medoid-based tree features for any sequence."""
    features = {'sequence_id': seq_id, 'label': label}
    
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    medoid_scores = {}
    for afp_type, info in medoid_info.items():
        stats = align_to_medoid(sequence, info['seq'], aligner)
        for k, v in stats.items():
            medoid_scores[f'medoid_{afp_type}_{k}'] = v
        features[f'medoid_{afp_type}_id'] = info['id']
    
    features.update(medoid_scores)
    
    # Top by identity
    m_idents = [(t.replace('medoid_', '').replace('_identity', ''), v)
                for t, v in medoid_scores.items() if t.endswith('_identity')]
    if m_idents:
        m_idents.sort(key=lambda x: x[1], reverse=True)
        features['medoid_top1_type'] = m_idents[0][0]
        features['medoid_top1_identity'] = m_idents[0][1]
        features['medoid_delta_top1_top2_ident'] = m_idents[0][1] - (m_idents[1][1] if len(m_idents) > 1 else 0.0)
    
    # Top by score
    m_scores = [(t.replace('medoid_', '').replace('_score', ''), v)
                for t, v in medoid_scores.items() if t.endswith('_score')]
    if m_scores:
        m_scores.sort(key=lambda x: x[1], reverse=True)
        features['medoid_top1_type_by_score'] = m_scores[0][0]
        features['medoid_delta_top1_top2_score'] = m_scores[0][1] - (m_scores[1][1] if len(m_scores) > 1 else 0.0)
    
    return features

def main():
    print("=== Extracting Medoid-Based Tree Features ===")
    print("WARNING: For CV, call get_medoids_for_types() with train-only data!")
    print("This standalone script uses all data (for demo/QC only).\n")
    
    afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']
    
    print("Finding medoids from UPGMA trees:")
    medoid_info = get_medoids_for_types(afp_types)
    
    if not medoid_info:
        print("Error: No medoids found!")
        return
    
    all_features = []
    
    # Process AFP sequences
    print("\nProcessing AFP sequences...")
    for record in SeqIO.parse("data/positives/afp_all_c90.faa", "fasta"):
        seq_id = record.id.split('|')[1] if '|' in record.id else record.id
        features = extract_medoid_features(str(record.seq), seq_id, medoid_info, 'AFP')
        all_features.append(features)
    print(f"  Processed {len([f for f in all_features if f['label'] == 'AFP'])} AFP sequences")
    
    # Process non-AFP sequences
    print("Processing non-AFP sequences...")
    sequences = list(SeqIO.parse("data/negatives/non_afp_c90.faa", "fasta"))
    for i, record in enumerate(sequences, 1):
        if i % 500 == 0:
            print(f"  Processed {i}/{len(sequences)}...")
        seq_id = record.id.split('|')[1] if '|' in record.id else record.id
        features = extract_medoid_features(str(record.seq), seq_id, medoid_info, 'NON_AFP')
        all_features.append(features)
    
    # Save
    df = pd.DataFrame(all_features)
    output_file = "data/features/medoid_tree_features.csv"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, index=False)
    
    print(f"\nSaved to {output_file}")
    print(f"Total: {len(df)} | AFP: {len(df[df['label']=='AFP'])} | Non-AFP: {len(df[df['label']=='NON_AFP'])}")
    print(f"Features per sequence: {len(df.columns)-2}")

if __name__ == "__main__":
    main()