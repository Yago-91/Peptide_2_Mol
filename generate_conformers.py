import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import argparse
import os

def process_molecule(data):
    smiles, mol_id, num_confs = data
    
    # RDKit expects strict string formatting
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    
    # Add hydrogens for accurate 3D geometry
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", str(mol_id))
    
    # Generate conformers using ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    
    # Optimize conformers with MMFF94 force field
    try:
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200)
    except Exception:
        pass # Skip optimization if force field fails for specific atoms
        
    return mol

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_excel", help="Excel file (.xlsx)")
    parser.add_argument("output_sdf", help="Output SDF file")
    parser.add_argument("--confs", type=int, default=50, help="Conformers per molecule")
    parser.add_argument("--cores", type=int, default=32, help="Number of CPU cores to use")
    args = parser.parse_args()

    print(f"Reading {args.input_excel}...")
    
    # Load Excel file. We drop completely empty rows.
    df = pd.read_excel(args.input_excel).dropna(how='all')
    
    tasks = []
    # Iterate through rows. iloc[:, 0] is column 1 (ID), iloc[:, 1] is column 2 (SMILES)
    for i in range(len(df)):
        try:
            mol_id = str(df.iloc[i, 0]).strip()
            smiles = str(df.iloc[i, 1]).strip()
            
            # Skip if either is "nan" (Pandas representation of empty cell)
            if mol_id.lower() != 'nan' and smiles.lower() != 'nan':
                tasks.append((smiles, mol_id, args.confs))
        except IndexError:
            continue

    print(f"Generating conformers using {args.cores} cores...")
    
    with multiprocessing.Pool(args.cores) as pool:
        results = pool.map(process_molecule, tasks)

    print("Writing to SDF...")
    writer = Chem.SDWriter(args.output_sdf)
    valid_count = 0
    for mol in results:
        if mol is not None:
            writer.write(mol)
            valid_count += 1
    writer.close()
    
    print(f"Successfully processed {valid_count} out of {len(tasks)} molecules.")

if __name__ == "__main__":
    main()