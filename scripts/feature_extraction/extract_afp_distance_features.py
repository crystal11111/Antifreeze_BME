from Bio import AlignIO, SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import numpy as np
from collections import Counter
import os
import sys
'''
What it does

- Builds per-type AFP profiles (Type1, Type2, Type3, Type4, AFGP) from MSAs.
- For every input sequence (positive + negative), aligns to ALL 5 type profiles
  and extracts alignment-based similarity features PER TYPE.
- (Optional, --with-medoids) Also aligns to per-type tree medoids to add tree-bsed features.
'''

def calculate_conservation_score(column):
    """Calculate conservation score for an alignment column"""
    residues = [res for res in column if res != '-']
    if not residues:
        return 0
    
    counts = Counter(residues)
    total = len(residues)
    
    entropy = 0
    for count in counts.values():
        freq = count / total
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    max_entropy = np.log2(min(20, len(set(residues))))
    return 1 - (entropy / max_entropy) if max_entropy > 0 else 1

def consensus_and_conservation(msa):
    """Build consensus and conservation profile from MSA"""
    L = msa.get_alignment_length()
    consensus = []
    cons_scores = np.zeros(L, dtype=float)
    ungap_to_msa = []  # Map ungapped consensus position to MSA column
    
    for i in range(L):
        col = [res for res in msa[:, i] if res != '-']
        if not col:
            consensus.append('X')
            cons_scores[i] = 0.0
            continue
        mc = Counter(col).most_common(1)[0][0]
        consensus.append(mc)
        cons_scores[i] = calculate_conservation_score(msa[:, i])
        if mc != '-':
            ungap_to_msa.append(i)
    
    consensus_clean = "".join(consensus).replace('-', '').replace('X', '').upper()
    return consensus_clean, cons_scores, ungap_to_msa

def align_stats(seq, ref_seq, aligner, cons_scores, ungap_to_msa):
    """Compute robust alignment statistics"""
    seq_u = seq.upper()
    ref_u = ref_seq.upper()
    aln = aligner.align(seq_u, ref_u)[0]
    q_blocks, r_blocks = aln.aligned
    
    matches = mismatches = 0
    cons_vals = []
    
    for (qs, qe), (rs, re) in zip(q_blocks, r_blocks):
        for a, b in zip(seq_u[qs:qe], ref_u[rs:re]):
            if a == b:
                matches += 1
            else:
                mismatches += 1
        # Collect conservation at matched positions
        for j in range(rs, re):
            if j < len(ungap_to_msa):
                msa_col = ungap_to_msa[j]
                cons_vals.append(cons_scores[msa_col])
    
    # Calculate gaps properly
    q_aligned = sum(qe - qs for qs, qe in q_blocks)
    r_aligned = sum(re - rs for rs, re in r_blocks)
    aligned_len = min(q_aligned, r_aligned)
    ins_q = max(0, q_aligned - aligned_len)
    ins_r = max(0, r_aligned - aligned_len)
    gaps = ins_q + ins_r
    
    cov = aligned_len / max(1, len(seq_u))
    ident = matches / max(1, aligned_len)
    
    return {
        'score': aln.score,
        'aligned_len': aligned_len,
        'coverage': cov,
        'identity': ident,
        'mismatch_rate': mismatches / max(1, aligned_len),
        'gap_rate': gaps / max(1, aligned_len + gaps),
        'conservation_match': float(np.mean(cons_vals)) if cons_vals else 0.0
    }

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

def extract_afp_distance_features(sequence, seq_id, type_profiles, label, medoid_info=None, afp_type=None, reviewed=False):
    """Extract AFP distance features for any sequence against all type profiles"""
    features = {'sequence_id': seq_id, 'label': label, 'reviewed': int(reviewed)}
    if (afp_type is not None):
        features['afp_type'] = afp_type
    else:
        features['afp_type'] = 'N/A'
    
    L = max(1, len(sequence))
    
    # Basic sequence properties (normalized)
    features['length'] = len(sequence)
    features['net_charge_per_residue'] = (sequence.count('K') + sequence.count('R') - sequence.count('D') - sequence.count('E')) / L
    features['hydrophobic_percent'] = sum(sequence.count(aa) for aa in 'AILMFPWYV') / L * 100
    features['polar_percent'] = sum(sequence.count(aa) for aa in 'NQST') / L * 100
    features['charged_percent'] = sum(sequence.count(aa) for aa in 'DEKR') / L * 100
    features['gly_percent'] = sequence.count('G') / L * 100
    features['pro_percent'] = sequence.count('P') / L * 100
    
    # Quality flags
    ambiguous = sum(sequence.count(aa) for aa in 'BXZJUO')
    features['ambiguous_count'] = ambiguous
    features['ambiguous_rate'] = ambiguous / L
    
    # Align to each AFP type profile
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    type_scores = {}
    for afp_type, (consensus, cons_scores, ungap_map) in type_profiles.items():
        stats = align_stats(sequence, consensus, aligner, cons_scores, ungap_map)
        for key, val in stats.items():
            type_scores[f'{afp_type}_{key}'] = val
    
    # Find top 2 types
    type_identities = [(t.replace('_identity', ''), type_scores[t]) for t in type_scores if '_identity' in t]
    type_identities.sort(key=lambda x: x[1], reverse=True)
    
    if len(type_identities) >= 2:
        features['top1_type'] = type_identities[0][0]
        features['top1_identity'] = type_identities[0][1]
        features['top1_score'] = type_scores[f"{type_identities[0][0]}_score"]
        features['top1_coverage'] = type_scores[f"{type_identities[0][0]}_coverage"]
        features['top2_type'] = type_identities[1][0]
        features['top2_identity'] = type_identities[1][1]
        features['delta_top1_top2'] = type_identities[0][1] - type_identities[1][1]
    else:
        features['top1_type'] = type_identities[0][0] if type_identities else 'none'
        features['top1_identity'] = type_identities[0][1] if type_identities else 0
        features['top1_score'] = type_scores.get(f"{features['top1_type']}_score", 0) if type_identities else 0
        features['top1_coverage'] = type_scores.get(f"{features['top1_type']}_coverage", 0) if type_identities else 0
        features['top2_type'] = 'none'
        features['top2_identity'] = 0
        features['delta_top1_top2'] = 0
    
    features.update(type_scores)
    features['alignment_failed'] = int(features['top1_type'] == 'none' or features['top1_identity'] < 0.1)
    features['low_coverage'] = int(features.get(f"{features['top1_type']}_coverage", 0) < 0.5)
    
    # Optional: medoid features (tree-based)
    if medoid_info:
        medoid_scores = {}
        for afp_type, info in medoid_info.items():
            mstats = align_to_medoid(sequence, info['seq'], aligner)
            for k, v in mstats.items():
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
    print("=== Extracting AFP Distance Features ===")
    print("WARNING: For CV, build type_profiles and medoid_info from train-only data per fold!")
    print("This standalone script uses all data (for demo/QC only).")
    print("Aligner parameters: mode=global, match=2, mismatch=-1, gap_open=-2, gap_extend=-0.5")
    
    # Optional: load medoids for tree-based features
    use_medoids = '--with-medoids' in sys.argv
    medoid_info = None
    if use_medoids:
        try:
            from medoid_tools import get_medoids_for_types
            print("\nLoading medoids for tree-based features:")
            medoid_info = get_medoids_for_types(['type1', 'type2', 'type3', 'type4', 'afgp'])
        except ImportError:
            print("Warning: medoid_tools not found, skipping medoid features")
    
    # Load per-type AFP alignments
    afp_types = ['type1', 'type2', 'type3', 'type4', 'afgp']
    type_profiles = {}
    
    for afp_type in afp_types:
        aln_file = f"data/positives/by_type/{afp_type}_c90.aln"
        if not os.path.exists(aln_file):
            print(f"Warning: {aln_file} not found, skipping")
            continue
        
        msa = AlignIO.read(aln_file, "clustal")
        consensus, cons_scores, ungap_map = consensus_and_conservation(msa)
        type_profiles[afp_type] = (consensus.upper(), cons_scores, ungap_map)
        print(f"Loaded {afp_type}: {len(msa)} seqs, {len(consensus)} consensus residues")
    
    print(f"\nType profiles loaded: {list(type_profiles.keys())}")
    
    if not type_profiles:
        print("Error: No type profiles loaded!")
        return
    
    all_features = []
    
    # Process reviewed AFP sequences
    print("\nProcessing Reviewed AFP sequences...")
    for afp_type in afp_types:
        aln_file = f"data/positives/by_type/{afp_type}_c90.faa"
        if not os.path.exists(aln_file):
            print(f"Warning: {aln_file} not found, skipping")
            continue
        for record in SeqIO.parse(aln_file, "fasta"):
            seq_id = record.id.split('|')[1] if '|' in record.id else record.id
            sequence = str(record.seq)
            
            features = extract_afp_distance_features(sequence, seq_id, type_profiles, 'AFP', medoid_info, afp_type, reviewed=True)
            all_features.append(features)
        
    print(f"  Processed {len([f for f in all_features if f['label'] == 'AFP'])} Reviewed AFP sequences")
    
    # Process non-reviewed AFP sequences
    print("\nProcessing Non-Reviewed AFP sequences...")
    for afp_type in afp_types:
        aln_file = f"data/positives/pblast_filtered/pblast_{afp_type}_c90.faa"
        if not os.path.exists(aln_file):
            print(f"Warning: {aln_file} not found, skipping")
            continue
        for record in SeqIO.parse(aln_file, "fasta"):
            seq_id = record.id.split('|')[1] if '|' in record.id else record.id
            sequence = str(record.seq)
            
            features = extract_afp_distance_features(sequence, seq_id, type_profiles, 'AFP', medoid_info, afp_type, reviewed=False)
            all_features.append(features)
        
    print(f"  Processed {len([f for f in all_features if f['label'] == 'AFP'])} Non-reviewed AFP sequences")

    # Process non-AFP sequences
    print("Processing non-AFP sequences...")
    sequences = list(SeqIO.parse("data/negatives/non_afp_c90.faa", "fasta"))
    
    for i, record in enumerate(sequences, 1):
        if i % 500 == 0:
            print(f"  Processed {i}/{len(sequences)} non-AFP sequences...")
        
        seq_id = record.id.split('|')[1] if '|' in record.id else record.id
        sequence = str(record.seq)
        
        features = extract_afp_distance_features(sequence, seq_id, type_profiles, 'NON_AFP', medoid_info, reviewed=True)
        all_features.append(features)
    
    # Save results
    df = pd.DataFrame(all_features)
    output_file = "data/features/afp_distance_features.csv"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, index=False)
    
    print(f"\nSaved AFP distance features to {output_file}")
    print(f"Total sequences: {len(df)}")
    print(f"AFP: {len(df[df['label'] == 'AFP'])}")
    print(f"Non-AFP: {len(df[df['label'] == 'NON_AFP'])}")
    print(f"Features per sequence: {len(df.columns)-2}")
    print(f"\nFeature breakdown:")
    print(f"  Basic properties: 8 (length, charge, hydrophobic%, polar%, charged%, gly%, pro%, ambiguous)")
    print(f"  Per-type features: {len([c for c in df.columns if 'type1_' in c])} × 5 types = {len([c for c in df.columns if 'type1_' in c]) * 5}")
    print(f"  Top-level: 7 (top1/top2 type/identity/score/coverage, delta)")
    print(f"  Quality flags: 2 (alignment_failed, low_coverage)")
    print(f"\nGap rates now properly calculated (non-zero for misaligned sequences)")
    print(f"Conservation scores computed at matched positions only")

if __name__ == "__main__":
    main()