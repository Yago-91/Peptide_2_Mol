import os
import argparse
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

# Mute RDKit warnings
RDLogger.DisableLog('rdApp.*')

# Set up command line arguments
parser = argparse.ArgumentParser(description="Extract hits from Excel and convert to PDBQT.")
parser.add_argument("-e", "--excel", required=True, help="Path to your Excel file with selected hits")
parser.add_argument("-s", "--sdf", required=True, help="Path to your aggregated master SDF file")
parser.add_argument("-o", "--outdir", default="Docking_Ready", help="Output directory for PDBQT files")
args = parser.parse_args()

# 1. Read the Excel file
print(f"Reading target IDs from {args.excel}...")
try:
    df = pd.read_excel(args.excel, header=0)
    target_ids = set(df.iloc[:, 1].dropna().astype(str).tolist())
except Exception as e:
    print(f"Error reading Excel file: {e}")
    exit(1)

print(f"Found {len(target_ids)} unique ZINC IDs to extract.")

# 2. Setup the output directory and the Modern Meeko Preparator
os.makedirs(args.outdir, exist_ok=True)
preparator = MoleculePreparation()

# 3. Scan the aggregated SDF
print(f"Scanning {args.sdf} and generating PDBQT files...")
suppl = Chem.SDMolSupplier(args.sdf)
processed_count = 0

for mol in suppl:
    if mol is None:
        continue
        
    zinc_id = mol.GetProp('_Name') if mol.HasProp('_Name') else None
    
    if zinc_id in target_ids:
        try:
            # Modern Meeko syntax (v0.6.0+)
            mol_setups = preparator.prepare(mol)
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
            
            if is_ok:
                out_path = os.path.join(args.outdir, f"{zinc_id}.pdbqt")
                with open(out_path, "w") as f:
                    f.write(pdbqt_string)
                processed_count += 1
                target_ids.remove(zinc_id) 
            else:
                print(f"Skipping {zinc_id}: {error_msg}")
                
        except Exception as e:
            print(f"Failed to prepare {zinc_id}: {e}")

print(f"\nSuccess! Generated {processed_count} ready-to-dock .pdbqt files in ./{args.outdir}/")

if len(target_ids) > 0:
    print(f"Warning: {len(target_ids)} IDs from your Excel were not found in the SDF.")