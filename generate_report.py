import csv
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit import RDLogger

parser = argparse.ArgumentParser(description="Extract, calculate, and deduplicate hits.")
parser.add_argument("-i", "--input", required=True, help="Raw SDF file")
parser.add_argument("-o", "--output", required=True, help="Output CSV report")
parser.add_argument("-s", "--sdf_out", required=True, help="Output SORTED SDF file")
args = parser.parse_args()

RDLogger.DisableLog('rdApp.*')

print(f"Extracting and deduplicating data from {args.input}...")

suppl = Chem.SDMolSupplier(args.input)
best_poses = {} 
score_tag_name = None

for mol in suppl:
    if mol is None:
        continue

    comp_id = mol.GetProp('_Name') if mol.HasProp('_Name') else "Unknown_ID"

    if score_tag_name is None:
        for prop in mol.GetPropNames():
            if "score" in prop.lower() or "fit" in prop.lower():
                score_tag_name = prop
                break
    
    if score_tag_name is None or not mol.HasProp(score_tag_name):
        continue

    score = float(mol.GetProp(score_tag_name))
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    mol_data = {
        "Compound_ID": comp_id,
        "Fit_Score": round(score, 4),
        "Mol_Weight": round(mw, 2),
        "LogP_Solubility": round(logp, 2),
        "Rotatable_Bonds": rot_bonds,
        "mol_obj": mol  # Keep the physical 3D object in memory
    }

    if comp_id not in best_poses:
        best_poses[comp_id] = mol_data
    else:
        if score > best_poses[comp_id]["Fit_Score"]:
            best_poses[comp_id] = mol_data

results = list(best_poses.values())
results.sort(key=lambda x: x['Fit_Score'], reverse=True)

print(f"Writing {len(results)} clean hits to {args.output} and {args.sdf_out}...")

# 1. Write the sorted CSV
with open(args.output, mode='w', newline='') as file:
    headers = ["Compound_ID", "Fit_Score", "Mol_Weight", "LogP_Solubility", "Rotatable_Bonds"]
    writer = csv.DictWriter(file, fieldnames=headers, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(results)

# 2. Write the perfectly synced, sorted SDF
writer_sdf = Chem.SDWriter(args.sdf_out)
for row in results:
    writer_sdf.write(row["mol_obj"])
writer_sdf.close()

print("Report and sorted SDF generation complete!")