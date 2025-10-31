- ncbi.fasta: NCBI datasets with query: "antifreeze protein"[Title/Abstract] AND (Actinopterygii OR Insecta)

- 1️⃣ Checking/cleaning your FASTA -> Use CD-HIT to remove redundant sequences (90% identity is standard).
    - cd-hit -i ncbi.fasta -o afp_nr90.fasta -c 0.9 -n 5 -d 0
- 2️⃣ Running MSA (Multiple Sequence Alignment)
    - mafft --localpair --maxiterate 1000 afp_nr90.fasta > afp_aligned.fasta
- 3️⃣ Detecting motifs (known and de novo)
    - run count_motifs.py
- 4️⃣ Saving features for later analysis (like regression modeling)
    - run combine_features.py