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
# Combine AFP sequences
cat data/raw/uniprotkb_fishAFP.fasta data/raw/uniprotkb_nonfishAFP.fasta > data/positives/afp_all_raw.faa

# Remove redundancy (90% identity clustering)
cd-hit -i data/positives/afp_all_raw.faa -o data/positives/afp_all_c90.faa -c 0.90 -n 5
```

### 2️⃣ Multiple Sequence Alignment
```bash
python conduct_MSA.py
```

### 3️⃣ Phylogenetic Analysis
```bash
python create_UPGMA.py
python plot_dnd.py
```

### 4️⃣ Motif Analysis
```bash
python count_motifs.py
```
### 5 training the HMM
```bash
python train_hmm.ipynb
```

## Scripts

- **conduct_MSA.py**: Performs multiple sequence alignment using ClustalW2
- **create_UPGMA.py**: Constructs UPGMA phylogenetic tree from alignment
- **plot_dnd.py**: Visualizes ClustalW guide tree from .dnd file
- **count_motifs.py**: Counts specific motifs (TxT, TxxT, TAA, TAP) in sequences
- **train_hmm.ipynb**: used to train the positive and negative HMM 
- **test_hmm.ipynb**: uses the models outputted by the HMM to predict 

## Requirements

- Python 3.x
- Biopython
- matplotlib
- ClustalW2 (must be in PATH)
- CD-HIT
- hmmlearn
- numpy