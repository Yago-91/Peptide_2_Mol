import argparse
import csv
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

# Setup command line arguments
parser = argparse.ArgumentParser(description="Cluster screening hits by Bemis-Murcko scaffolds.")
parser.add_argument("-i", "--input", required=True, help="Your Sorted SDF file from the reporter script.")
parser.add_argument("-o", "--output", required=True, help="Output CSV for the scaffold analysis.")
args = parser.parse_args()

RDLogger.DisableLog('rdApp.*')

print(f"Analyzing chemical scaffolds in {args.input}...")

suppl = Chem.SDMolSupplier(args.input)
scaffold_data = {}

for mol in suppl:
    if mol is None:
        continue
    
    comp_id = mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown_ID"
    
    # Try to extract the score (assuming it was saved by your reporter script)
    score = 0.0
    for prop in mol.GetPropNames():
        if "score" in prop.lower() or "fit" in prop.lower():
            score = float(mol.GetProp(prop))
            break

    # 1. Strip the side-chains to isolate the core Bemis-Murcko Scaffold
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        # Convert the 3D core into a 1D SMILES string so we can group identical ones
        core_smiles = Chem.MolToSmiles(core)
    except Exception:
        continue # Skip if RDKit cannot extract a scaffold

    # 2. Group and count the scaffolds
    if core_smiles not in scaffold_data:
        scaffold_data[core_smiles] = {
            "Frequency": 1,
            "Best_Score": score,
            "Best_Representative_ID": comp_id
        }
    else:
        scaffold_data[core_smiles]["Frequency"] += 1
        # Update the best score and ID if this specific molecule scored higher
        if score > scaffold_data[core_smiles]["Best_Score"]:
            scaffold_data[core_smiles]["Best_Score"] = score
            scaffold_data[core_smiles]["Best_Representative_ID"] = comp_id

# 3. Sort the scaffolds by how frequently they appear in your hits
sorted_scaffolds = sorted(scaffold_data.items(), key=lambda x: x[1]['Frequency'], reverse=True)

print(f"Distilled hits into {len(sorted_scaffolds)} unique chemical cores.")
print(f"Writing cluster report to {args.output}...")

with open(args.output, mode='w', newline='') as file:
    headers = ["Scaffold_SMILES", "Frequency", "Best_Score", "Best_Representative_ID"]
    writer = csv.DictWriter(file, fieldnames=headers)
    writer.writeheader()
    
    for smiles, data in sorted_scaffolds:
        writer.writerow({
            "Scaffold_SMILES": smiles,
            "Frequency": data["Frequency"],
            "Best_Score": round(data["Best_Score"], 4),
            "Best_Representative_ID": data["Best_Representative_ID"]
        })

print("Clustering complete! You can now review your most dominant fragments.")