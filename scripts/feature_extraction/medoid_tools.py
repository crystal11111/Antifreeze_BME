from Bio import Phylo, SeqIO

def _clean_id(name: str) -> str:
    return name.split('|')[1] if '|' in name else name

def medoid_from_newick(tree_file: str) -> str:
    """Return the sequence ID of the medoid: node with min total path length to all others."""
    tree = Phylo.read(tree_file, "newick")
    
    # Check branch lengths
    for clade in tree.find_clades():
        if clade.branch_length is None:
            raise ValueError(f"{tree_file} has missing branch lengths")
    
    terms = tree.get_terminals()
    ids = [_clean_id(t.name) for t in terms]
    
    totals = []
    for i, ti in enumerate(terms):
        s = sum(tree.distance(ti, tj) for j, tj in enumerate(terms) if i != j)
        totals.append(s)
    
    medoid_idx = min(range(len(terms)), key=lambda k: totals[k])
    return ids[medoid_idx]

def load_medoid_sequence(fasta_path: str, medoid_id: str) -> str:
    """Return amino-acid sequence string for a medoid ID from a FASTA."""
    for rec in SeqIO.parse(fasta_path, "fasta"):
        rid = _clean_id(rec.id)
        if rid == medoid_id:
            return str(rec.seq).upper()
    raise ValueError(f"Medoid id {medoid_id} not found in {fasta_path}")

def get_medoids_for_types(afp_types, tree_fmt="data/positives/by_type/{type}_c90.nwk", 
                          fasta_fmt="data/positives/by_type/{type}_c90.faa"):
    """Get medoid info (ID + sequence) for all AFP types. Call with train-only data in CV."""
    medoid_info = {}
    for t in afp_types:
        tree_file = tree_fmt.format(type=t)
        fasta_file = fasta_fmt.format(type=t)
        try:
            medoid_id = medoid_from_newick(tree_file)
            medoid_seq = load_medoid_sequence(fasta_file, medoid_id)
            medoid_info[t] = {"id": medoid_id, "seq": medoid_seq}
            print(f"  {t}: medoid = {medoid_id}")
        except (FileNotFoundError, ValueError) as e:
            print(f"  {t}: skipped ({e})")
    return medoid_info
