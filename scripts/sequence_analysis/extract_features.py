# Potential features we can use to analyze the sequence are as follows:
# 1) frequency of identifed motifs - from MSA
# 2) Sequence length (as Type 3 proteins tend to be smaller)
# 3) Hydrophobicity
# 4) UPGMA distance from nearest AFP clade
# 5) HMM negative and positive log probabilities
# 6) PSSM from PSI-BLAST
# 7) Steric hindrance values

import numpy as np

# Taken from the qiagen bioninformatics manual at:
# https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/2200/index.php?manual=BE_Protein_hydrophobicity.html
hydropathy = {
    'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,
    'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,
    'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3
}

def hydrophobicity_score(seq):
    return np.mean([hydropathy[a] for a in seq if a in hydropathy])


# Taken from the proteins and proteomics organization headed by Professor Simpson
# https://proteinsandproteomics.org/content/free/tables_1/table08.pdf
steric_volume = {
    'A': 67, 'C': 86, 'D': 91, 'E': 109, 'F': 135, 'G': 48, 'H': 118,
    'I': 124,'K': 135,'L': 124,'M': 124,'N': 96,'P': 90,'Q': 114,
    'R': 148,'S': 73,'T': 93,'V': 105,'W': 163,'Y': 141
}

def steric_hindrance_score(seq):
    return sum(steric_volume[a] for a in seq if a in steric_volume) / len(seq)

