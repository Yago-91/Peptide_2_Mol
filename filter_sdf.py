import pandas as pd
from rdkit import Chem

# --- 1. SET YOUR FILE NAMES HERE ---
excel_file = 'Filtered_hits.xlsx'
input_sdf = 'All_unfiltered_Hits.sdf'

output_sdf = 'All_filtered_hits.sdf'
output_smiles = 'All_filtered_hits_smiles.txt'

# --- 2. EXTRACT IDs FROM EXCEL ---
print("Reading Excel file...")
# Pandas automatically uses the first row as a header.
df = pd.read_excel(excel_file)

# .iloc[:, 1] grabs the second column (index 1). 
# We convert it to a 'set' because looking up items in a set is lightning fast in Python.
hit_ids = set(df.iloc[:, 1].astype(str).tolist())
print(f"Loaded {len(hit_ids)} unique target IDs from Excel.")

# --- 3. FILTER THE SDF ---
print("Scanning SDF file...")
supplier = Chem.SDMolSupplier(input_sdf)
sdf_writer = Chem.SDWriter(output_sdf)
smiles_writer = open(output_smiles, 'w')

match_count = 0

for mol in supplier:
    if mol is not None:
        # Extract the molecule's ID. In standard SDFs, this is the title block ('_Name')
        if mol.HasProp('_Name'):
            mol_id = str(mol.GetProp('_Name')).strip()
        else:
            mol_id = "UNKNOWN"
        
        # If the ID from the SDF matches an ID in our Excel set, keep it
        if mol_id in hit_ids:
            # 1. Write to the new filtered SDF
            sdf_writer.write(mol)
            
            # 2. Generate and write the SMILES
            smiles = Chem.MolToSmiles(mol)
            smiles_writer.write(f"{mol_id}\t{smiles}\n")
            
            match_count += 1

# --- 4. CLEANUP ---
sdf_writer.close()
smiles_writer.close()

print(f"Done! Successfully extracted {match_count} compounds.")
print(f"Saved to '{output_sdf}' and '{output_smiles}'.")