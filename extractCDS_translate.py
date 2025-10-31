from Bio import Entrez, SeqIO
from Bio.Seq import Seq

Entrez.email = "crystal031230@gmail.com"

term = '(antifreeze[Title/Abstract]) AND (Actinopterygii[Organism] OR Insecta[Organism]) AND (mrna[Filter])'
handle = Entrez.esearch(db="nucleotide", term=term, retmax=200)
ids = Entrez.read(handle)['IdList']

fasta_out = open("mrna_raw.fasta", "w")
for gi in ids:
    rec = SeqIO.read(Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="text"), "genbank")
    # If CDS is annotated, translate those
    cds_features = [f for f in rec.features if f.type == "CDS"]
    for i, f in enumerate(cds_features):
        seq = f.extract(rec.seq)
        prot = seq.translate(to_stop=True)
        header = f">{rec.id}|cds{i}|{f.qualifiers.get('product',['NA'])[0]}"
        fasta_out.write(f"{header}\n{str(prot)}\n")
fasta_out.close()
