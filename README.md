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
# Extract MSA-based features (conservation, consensus identity)
python scripts/feature_extraction/extract_msa_features.py

# Extract phylogenetic tree features (branch lengths, distances)
python scripts/feature_extraction/extract_tree_features.py

# Basic sequence properties (hydrophobicity/steric hindrance)
python scripts/sequence_analysis/extract_features.py
```


## Scripts

### Data Processing
- **split_afp_by_type.py**: Splits fish AFPs into separate files by type (Type1, Type2, Type3, Type4, AFGP)
- **process_afp_types.py**: Processes each AFP type through cd-hit clustering and MSA

### Sequence Analysis
- **conduct_MSA.py**: Performs multiple sequence alignment using ClustalW
- **create_UPGMA.py**: Constructs UPGMA phylogenetic tree from alignment
- **create_UPGMA_by_type.py**: Generates UPGMA trees for all AFP types
- **plot_dnd.py**: Visualizes ClustalW guide tree from .dnd file
- **count_motifs.py**: Counts specific motifs (TxT, TxxT, TAA, TAP) in sequences
- **extract_features.py**: Calculates sequence hydrophobicity and steric hindrance
- **common_motifs.py**: Gets common motifs from an .aln file

### Feature Extraction
- **extract_msa_features.py**: Extracts MSA-based features (conservation scores, consensus identity)
- **extract_tree_features.py**: Extracts phylogenetic tree features (branch lengths, distances) 

### Machine Learning
- **train_hmm.ipynb**: Trains positive and negative HMM models
- **test_hmm.ipynb**: Uses trained models for AFP prediction

## Requirements

- Python 3.x
- Biopython
- matplotlib
- pandas
- numpy
- hmmlearn
- ClustalW (must be in PATH)
- CD-HIT
- Jupyter Notebook