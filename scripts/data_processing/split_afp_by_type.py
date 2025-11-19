import pandas as pd
from Bio import SeqIO
import os
from collections import Counter

# Read metadata
meta_df = pd.read_csv("data/meta/uniprotkb_fishAFP.tsv", sep="\t")

# Create type mapping from protein names
type_mapping = {}
for _, row in meta_df.iterrows():
    entry_id = row['Entry']
    protein_name = row['Protein names']
    
    # Determine AFP type from protein name
    if 'Type-2' in protein_name or 'Type II' in protein_name:
        afp_type = 'Type2'
    elif 'Type-3' in protein_name or 'Type III' in protein_name:
        afp_type = 'Type3'
    elif 'Type-4' in protein_name or 'Type IV' in protein_name:
        afp_type = 'Type4'
    elif 'glycoprotein' in protein_name or 'AFGP' in protein_name:
        afp_type = 'AFGP'
    else:
        afp_type = 'Type1'  # Default for others
    
    type_mapping[entry_id] = afp_type

# Define typical organisms and length ranges for each type
type_criteria = {
    'Type1': {'organisms': ['Pseudopleuronectes americanus', 'Myoxocephalus'], 'length_range': (30, 45)},  # 3.3-4.5 kD
    'Type2': {'organisms': ['Hemitripterus americanus', 'Osmerus mordax'], 'length_range': (120, 140)},
    'Type3': {'organisms': ['Zoarces americanus', 'Anarhichas lupus'], 'length_range': (50, 65)},  # ~6 kD
    'Type4': {'organisms': ['Myoxocephalus octodecemspinosus', 'Gadus morhua'], 'length_range': (100, 120)},  # ~12 kD
    'AFGP': {'organisms': ['Eleginus gracilis', 'Notothenia', 'Pagothenia'], 'length_range': (20, 35)}  # 2.6-3.3 kD
}

# Apply reclassifications based on organism and length criteria
print("--- Applying Automatic Reclassifications --")
reclassifications = {}

for _, row in meta_df.iterrows():
    entry_id = row['Entry']
    current_type = type_mapping[entry_id]
    protein_name = row['Protein names']
    organism = row['Organism']
    length = row['Length']
    
    # Skip if has explicit type naming
    has_explicit_name = (
        (current_type == 'Type1' and ('Type-1' in protein_name or 'Type 1' in protein_name or 'Type I' in protein_name)) or
        (current_type == 'Type2' and ('Type-2' in protein_name or 'Type 2' in protein_name or 'Type II' in protein_name)) or
        (current_type == 'Type3' and ('Type-3' in protein_name or 'Type 3' in protein_name or 'Type III' in protein_name)) or
        (current_type == 'Type4' and ('Type-4' in protein_name or 'Type 4' in protein_name or 'Type IV' in protein_name)) or
        (current_type == 'AFGP' and ('glycoprotein' in protein_name or 'AFGP' in protein_name))
    )
    
    if has_explicit_name:
        continue
    
    # Check for better type match
    best_alternative = None
    best_score = 0
    
    for alt_type in ['Type1', 'Type2', 'Type3', 'Type4', 'AFGP']:
        if alt_type == current_type:
            continue
            
        alt_score = 0
        
        # Check organism match (70%)
        for typical_org in type_criteria[alt_type]['organisms']:
            if typical_org in organism:
                alt_score = 70
                break
        
        # Check length range (+15% or 30% alone)
        alt_length_range = type_criteria[alt_type]['length_range']
        if alt_length_range[0] <= length <= alt_length_range[1]:
            if alt_score > 0:
                alt_score += 15
            else:
                alt_score = 30
        
        if alt_score > best_score:
            best_score = alt_score
            best_alternative = alt_type
    
    # Apply reclassification if score >= 80%
    if best_alternative and best_score >= 80:
        reclassifications[entry_id] = (current_type, best_alternative, best_score)
        type_mapping[entry_id] = best_alternative

print(f"Reclassified {len(reclassifications)} proteins:")
for entry_id, (old_type, new_type, score) in reclassifications.items():
    print(f"  {entry_id}: {old_type} → {new_type} ({score}%)")

# Create output directory
os.makedirs("data/positives/by_type", exist_ok=True)

# Remove old files
for afp_type in ['type1', 'type2', 'type3', 'type4', 'afgp']:
    filepath = f"data/positives/by_type/{afp_type}_fish.fasta"
    if os.path.exists(filepath):
        os.remove(filepath)

# Initialize file handles with corrected types
type_files = {}
for afp_type in set(type_mapping.values()):
    type_files[afp_type] = open(f"data/positives/by_type/{afp_type.lower()}_fish.fasta", "w")

# Parse FASTA and split by corrected type
for record in SeqIO.parse("data/raw/uniprotkb_fishAFP.fasta", "fasta"):
    entry_id = record.id.split("|")[1]
    afp_type = type_mapping.get(entry_id, 'Type1')
    SeqIO.write(record, type_files[afp_type], "fasta")

# Close files
for f in type_files.values():
    f.close()

# Print final distribution
print("\n--- Final AFP Type Distribution ---")
final_counts = Counter(type_mapping.values())
total = len(type_mapping)

for afp_type in ['Type1', 'Type2', 'Type3', 'Type4', 'AFGP']:
    if afp_type in final_counts:
        count = final_counts[afp_type]
        percentage = (count / total) * 100
        print(f"{afp_type}: {count} sequences ({percentage:.1f}%)")

print(f"\nTotal: {total} sequences processed")