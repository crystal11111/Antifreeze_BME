# Antifreeze Protein Analysis Pipeline

A bioinformatics pipeline for analyzing antifreeze proteins (AFPs) using multiple sequence alignment, phylogenetic analysis, and machine learning.

## Project Structure

```
Antifreeze_BME/
├── data/
│   ├── features/      # Extracted features for ML
│   ├── motifs/        # Common motifs from our data
│   ├── meta/          # UniProt metadata files
│   ├── raw/           # Original FASTA files
│   ├── positives/     # AFP sequences
│   │   └── by_type/   # Type-specific AFP files
│   └── negatives/     # Non-AFP sequences
├── scripts/
│   ├── data_processing/     # Data preparation scripts
│   ├── sequence_analysis/   # MSA, phylogeny, motifs
│   ├── feature_extraction/  # Feature extraction tools
│   └── machine_learning/    # HMM training/testing
├── model_results/       # HMM model outputs
└── README.md
```

## Data Sources

### UniProt Datasets
- **uniprotkb_nonfishAFP.tsv**: Non-fish antifreeze proteins
  - Query: `reviewed:true AND keyword:KW-0047 AND NOT taxonomy_id:7898`
- **uniprotkb_fishAFP.tsv**: Fish antifreeze proteins  
  - Query: `reviewed:true AND keyword:KW-0047 AND taxonomy_id:7898`
- **uniprotkb_allfishNonAFP.tsv**: Fish proteins (non-AFP)
  - Query: `reviewed:true AND NOT keyword:KW-0047 AND taxonomy_id:7898`
- **uniprotkb_allAFP.tsv**: All antifreeze proteins
  - Query: `reviewed:true AND keyword:KW-0047`

## Workflow

### 1️⃣ Data Preparation
```bash
# Split fish AFPs by type (Type1, Type2, Type3, Type4, AFGP)
python scripts/data_processing/split_afp_by_type.py

# Process each AFP type (cd-hit clustering + MSA)
python scripts/data_processing/process_afp_types.py
```

### 2️⃣ Phylogenetic Analysis
```bash
# Generate UPGMA trees for all AFP types
python scripts/sequence_analysis/create_UPGMA_by_type.py

# For individual types (modify input file path as needed):
python scripts/sequence_analysis/conduct_MSA.py
python scripts/sequence_analysis/create_UPGMA.py
python scripts/sequence_analysis/plot_dnd.py
```

### 3️⃣ Motif Analysis
```bash
# Counts motifs found in research papers
python scripts/sequence_analysis/count_motifs.py

# Extracts common motifs from our sequences
python scripts/sequence_analysis/common_motifs.py
```

### 4️⃣ Machine Learning
```bash
# Train HMM models
jupyter notebook scripts/machine_learning/train_hmm.ipynb

# Test HMM predictions
jupyter notebook scripts/machine_learning/test_hmm.ipynb
```

### 5️⃣ Feature Extraction
```bash
# Main classifier features (for AFP detection)
python scripts/feature_extraction/extract_afp_distance_features.py
python scripts/feature_extraction/extract_tree_features.py
```

## Feature Datasets

### For AFP Detection (Classification)
**Primary features** - Use these for training AFP vs non-AFP classifiers:

1. **afp_distance_features.csv** (4,738 sequences × 53 features)
   - Aligns all sequences against AFP type profiles
   - 8 basic properties: length, charge, amino acid composition
   - 35 per-type features: alignment score, identity, coverage, gap rate, conservation match (7 metrics × 5 types)
   - 7 discriminators: top1/top2 type, identity, score, coverage, delta
   - 3 quality flags: alignment_failed, low_coverage, ambiguous_rate

2. **medoid_tree_features.csv** (4,738 sequences × ~25 features)
   - Aligns all sequences to medoid (representative) of each AFP type
   - Per-type: score, identity, coverage, gap rate
   - Top medoid matches and deltas

## Scripts

### Data Processing
- **split_afp_by_type.py**: Splits fish AFPs by type (Type1, Type2, Type3, Type4, AFGP)
- **process_afp_types.py**: cd-hit clustering + MSA for each AFP type

### Sequence Analysis
- **conduct_MSA.py**: Performs multiple sequence alignment using ClustalW
- **create_UPGMA.py**: Constructs UPGMA phylogenetic tree from alignment
- **create_UPGMA_by_type.py**: Generates UPGMA trees for all AFP types
- **plot_dnd.py**: Visualizes ClustalW guide tree from .dnd file
- **count_motifs.py**: Counts specific motifs (TxT, TxxT, TAA, TAP) in sequences
- **common_motifs.py**: Extracts common motifs from alignments

### Feature Extraction
- **extract_afp_distance_features.py**: Per-type alignment features for all sequences
- **extract_tree_features.py**: Medoid-based phylogenetic distance features
- **medoid_tools.py**: Identifies medoid sequences from UPGMA trees
- **extract_msa_features.py**: Within-type conservation features (AFP-only)

### Machine Learning
- **train_hmm.ipynb**: Trains positive and negative HMM models
- **test_hmm.ipynb**: Uses trained models for AFP prediction

## Requirements

- Python 3.x
- Biopython
- matplotlib, pandas, numpy
- hmmlearn
- ClustalW (must be in PATH)
- CD-HIT
- Jupyter Notebook

## Quick Start for ML

```bash
# 1. Extract features
python scripts/feature_extraction/extract_afp_distance_features.py
python scripts/feature_extraction/extract_tree_features.py

# 2. Load features for training
import pandas as pd
df_dist = pd.read_csv('data/features/afp_distance_features.csv')
df_tree = pd.read_csv('data/features/medoid_tree_features.csv')

# 3. Merge on sequence_id
df = df_dist.merge(df_tree, on=['sequence_id', 'label'])

# 4. Train your classifier
X = df.drop(['sequence_id', 'label'], axis=1)
y = (df['label'] == 'AFP').astype(int)
```