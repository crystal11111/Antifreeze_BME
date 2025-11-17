import re
from Bio import SeqIO

motifs = {
    "TxT": re.compile(r"T[A-Z]T"),
    "TxxT": re.compile(r"T[A-Z]{2}T"),
    "TAA": re.compile(r"TAA"),
    "TAP": re.compile(r"TAP"),
}

print("accession\tlength\tTxT\tTxxT\tTAA\tTAP")
for rec in SeqIO.parse("afp_trimmed.fasta", "fasta"):
    s = str(rec.seq)
    counts = [len(m.findall(s)) for m in motifs.values()]
    print(f"{rec.id}\t{len(s)}\t" + "\t".join(map(str, counts)))
