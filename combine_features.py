import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

records = []
for rec in SeqIO.parse("afp_trimmed.fasta", "fasta"):
    seq = str(rec.seq)
    analyzed = ProteinAnalysis(seq)
    features = {
        "accession": rec.id,
        "length": len(seq),
        "percent_Thr": seq.count("T") / len(seq),
        "percent_Ala": seq.count("A") / len(seq),
        "hydrophobicity": analyzed.gravy(),
        "charge": analyzed.charge_at_pH(7.0),
        "aromaticity": analyzed.aromaticity(),
    }
    records.append(features)

pd.DataFrame(records).to_csv("sequence_features.csv", index=False)
