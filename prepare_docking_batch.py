import os
import argparse
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from meeko import MoleculePreparation
from meeko import PDBQTMolecule

# Mute RDKit warnings for messy un-filtered hits
RDLogger.DisableLog('rdApp.*')

# Set up command line arguments
parser = argparse.ArgumentParser(description="Extract hits from Excel and convert to PDBQT.")
parser.add_argument("-e", "--excel", required=True, help="Path to your Excel file with selected hits")
parser.add_argument("-s", "--sdf", required=True, help="Path to your aggregated master SDF file")
parser.add_argument("-o", "--outdir", default="Docking_Ready", help="Output directory for PDBQT files")
args = parser.parse_args()

# 1. Read the Excel file to get the target ZINC IDs
print(f"Reading target IDs from {args.excel}...")
# Assuming Column 1 is Campaign and Column 2 is ZINC ID based on your description
try:
    df = pd.read_excel(args.excel, header=True)
    # Extract the second column (index 1) and convert to a set for lightning-fast lookups
    target_ids = set(df.iloc[:, 1].dropna().astype(str).tolist())
except Exception as e:
    print(f"Error reading Excel file: {e}")
    exit(1)

print(f"Found {len(target_ids)} unique ZINC IDs to extract.")

# 2. Setup the output directory and the Meeko Preparator
os.makedirs(args.outdir, exist_ok=True)
preparator = MoleculePreparation()

# 3. Scan the aggregated SDF and convert matches directly to PDBQT
print(f"Scanning {args.sdf} and generating PDBQT files...")
suppl = Chem.SDMolSupplier(args.sdf)
processed_count = 0

for mol in suppl:
    if mol is None:
        continue
        
    zinc_id = mol.GetProp('_Name') if mol.HasProp('_Name') else None
    
    if zinc_id in target_ids:
        try:
            # Prepare the molecule for AutoDock (assigns charges, finds rotatable bonds)
            preparator.prepare(mol)
            pdbqt_string = preparator.write_pdbqt_string()
            
            # Write directly to the final file
            out_path = os.path.join(args.outdir, f"{zinc_id}.pdbqt")
            with open(out_path, "w") as f:
                f.write(pdbqt_string)
                
            processed_count += 1
            # Remove the ID from the set so we don't process duplicates if they exist
            target_ids.remove(zinc_id) 
            
        except Exception as e:
            print(f"Failed to prepare {zinc_id}: {e}")

print(f"\nSuccess! Generated {processed_count} ready-to-dock .pdbqt files in ./{args.outdir}/")

if len(target_ids) > 0:
    print(f"Warning: {len(target_ids)} IDs from your Excel were not found in the SDF.")