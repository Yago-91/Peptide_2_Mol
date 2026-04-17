import pandas as pd
from rdkit import Chem

# --- 1. SET YOUR FILE NAMES HERE ---
excel_file = 'Filtered_hits.xlsx'
input_sdf = 'All_Unfiltered_Hits.sdf'

output_sdf = 'All_filtered_hits.sdf'
output_smiles = 'All_filtered_hits_smiles.txt'

# --- 2. EXTRACT IDs FROM EXCEL ---
print("Reading Excel file...")
df = pd.read_excel(excel_file)
hit_ids = set(df.iloc[:, 1].astype(str).tolist())
print(f"Loaded {len(hit_ids)} unique target IDs from Excel.")

# --- 3. FILTER THE SDF (WITH SANITIZATION FIX) ---
print("Scanning SDF file...")

# Load with sanitize=False so RDKit doesn't immediately reject the 5-valent Nitrogens
supplier = Chem.SDMolSupplier(input_sdf, sanitize=False)
sdf_writer = Chem.SDWriter(output_sdf)
smiles_writer = open(output_smiles, 'w')

match_count = 0

for mol in supplier:
    if mol is not None:
        try:
            # Perform partial sanitization: Update properties loosely, then sanitize 
            # everything EXCEPT the strict valence/property rules.
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
            
            # Extract the molecule's ID
            if mol.HasProp('_Name'):
                mol_id = str(mol.GetProp('_Name')).strip()
            else:
                mol_id = "UNKNOWN"
            
            # If it's a hit, write it
            if mol_id in hit_ids:
                sdf_writer.write(mol)
                smiles = Chem.MolToSmiles(mol)
                smiles_writer.write(f"{mol_id}\t{smiles}\n")
                match_count += 1
                
        except Exception as e:
            # This catches any truly broken molecules that still fail
            print(f"Skipped a severely malformed molecule. Error: {e}")

# --- 4. CLEANUP ---
sdf_writer.close()
smiles_writer.close()

print(f"Done! Successfully extracted {match_count} compounds.")
print(f"Saved to '{output_sdf}' and '{output_smiles}'.")