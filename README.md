- uniprotkb_nonfishAFP.tsv
- reviewed:true AND keyword:KW-0047 AND NOT taxonomy_id:7898

- uniprotkb_fishAFP.tsv
- reviewed:true AND keyword:KW-0047 AND taxonomy_id:7898

- uniprotkb_allfishNonAFP.tsv
- reviewed:true AND NOT keyword:KW-0047 AND taxonomy_id:7898

- uniprotkb_allAFP.tsv
- reviewed:true AND keyword:KW-0047

cat data/raw/uniprotkb_fishAFP.fasta data/raw/uniprotkb_nonfishAFP.fasta > data/positives/afp_all_raw.faa
cd-hit -i data/positives/afp_all_raw.faa -o data/positives/afp_all_c90.faa -c 0.90 -n 5

- 2️⃣ Running MSA (Multiple Sequence Alignment)
"""
cat data/raw/uniprotkb_fishAFP.fasta data/raw/uniprotkb_nonfishAFP.fasta > data/positives/afp_all_raw.faa
cd-hit -i data/positives/afp_all_raw.faa -o data/positives/afp_all_c90.faa -c 0.90 -n 5"""