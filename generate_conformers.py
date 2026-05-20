import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import argparse
import sys

def process_molecule(data):
    smiles, mol_id, num_confs = data
    
    try:
        # 1. Parse the SMILES
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
            
        # 2. Strict Organic Filter (CHNOPS + B, Si, Halogens)
        # Atomic numbers for: H, B, C, N, O, F, Si, P, S, Cl, Br, I
        allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in allowed_atomic_nums:
                # If it has Pt, Sr, As, etc., silently drop it
                return None 
                
        # 3. Add hydrogens for accurate 3D geometry
        mol = Chem.AddHs(mol)
        mol.SetProp("_Name", str(mol_id))
        
        # 4. Generate conformers using ETKDG
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.clearConfs = True
        
        # Attempt to embed. If it returns -1, it failed to generate 3D coordinates.
        if AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params) == -1:
            return None
        
        # 5. Optimize conformers with MMFF94 force field
        # We loop through manually so if one pose fails, we keep the rest
        for conf in mol.GetConformers():
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=conf.GetId(), maxIters=200)
            except Exception:
                continue # Skip optimization for this specific pose if force field fails
            
        return mol

    except Exception:
        # If anything triggers a severe C++ crash or Python exception, drop the molecule
        # and prevent the multiprocessing thread from hanging.
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_excel", help="Excel file (.xlsx)")
    parser.add_argument("output_sdf", help="Output SDF file")
    parser.add_argument("--confs", type=int, default=50, help="Conformers per molecule")
    parser.add_argument("--cores", type=int, default=32, help="Number of CPU cores to use")
    args = parser.parse_args()

    print(f"Reading {args.input_excel}...")
    
    df = pd.read_excel(args.input_excel).dropna(how='all')
    
    tasks = []
    for i in range(len(df)):
        try:
            mol_id = str(df.iloc[i, 0]).strip()
            smiles = str(df.iloc[i, 1]).strip()
            
            if mol_id.lower() != 'nan' and smiles.lower() != 'nan':
                tasks.append((smiles, mol_id, args.confs))
        except IndexError:
            continue

    print(f"Generating conformers for {len(tasks)} valid inputs using {args.cores} cores...")
    print("Molecules containing heavy metals/unsupported elements will be safely ignored.")
    
    with multiprocessing.Pool(args.cores) as pool:
        results = pool.map(process_molecule, tasks)

    print("Writing valid conformers to SDF...")
    writer = Chem.SDWriter(args.output_sdf)
    valid_mols = 0
    total_poses = 0
    
    for mol in results:
        if mol is not None and mol.GetNumConformers() > 0:
            for conf in mol.GetConformers():
                writer.write(mol, confId=conf.GetId())
                total_poses += 1
            valid_mols += 1
            
    writer.close()
    
    print(f"Successfully processed {valid_mols} molecules, generating {total_poses} total 3D poses.")

if __name__ == "__main__":
    main()