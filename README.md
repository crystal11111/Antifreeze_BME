# Antifreeze Protein Analysis Pipeline

A bioinformatics pipeline for analyzing antifreeze proteins (AFPs) using multiple sequence alignment and phylogenetic analysis.

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
python split_afp_by_type.py

# Process each AFP type (cd-hit clustering + MSA)
python process_afp_types.py
```

### 2️⃣ Individual Type Analysis
```bash
# For specific AFP types, use existing scripts:
python conduct_MSA.py    # (modify input file path as needed)
python create_UPGMA.py   # (modify input file path as needed)
python plot_dnd.py       # (modify input file path as needed)
```

### 3️⃣ Motif Analysis
```bash
python count_motifs.py
```
### 5 training the HMM
```bash
python train_hmm.ipynb
```

## Scripts

- **split_afp_by_type.py**: Splits fish AFPs into separate files by type (Type1, Type2, Type3, Type4, AFGP)
- **process_afp_types.py**: Processes each AFP type through cd-hit clustering and MSA
- **conduct_MSA.py**: Performs multiple sequence alignment using ClustalW
- **create_UPGMA.py**: Constructs UPGMA phylogenetic tree from alignment
- **plot_dnd.py**: Visualizes ClustalW guide tree from .dnd file
- **count_motifs.py**: Counts specific motifs (TxT, TxxT, TAA, TAP) in sequences
- **train_hmm.ipynb**: used to train the positive and negative HMM 
- **test_hmm.ipynb**: uses the models outputted by the HMM to predict 

## Requirements

- Python 3.x
- Biopython
- matplotlib
- pandas
- ClustalW (must be in PATH)
- CD-HIT
- hmmlearn
- numpy