import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import argparse
import signal
import sys

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Molecule C++ embedding timed out")

def process_molecule(data):
    smiles, mol_id, num_confs = data
    
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(60)
    
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
            
        allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in allowed_atomic_nums:
                return None 
                
        if mol.GetNumHeavyAtoms() > 60:
            return None
                
        mol = Chem.AddHs(mol)
        
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.clearConfs = True
        params.maxIterations = 1000 
        
        if AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params) == -1:
            return None
        
        for conf in mol.GetConformers():
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=conf.GetId(), maxIters=200)
            except Exception:
                continue 
        
        # RETURN BOTH the molecule and the ID safely to avoid pickling loss
        return (mol, str(mol_id).strip())

    except TimeoutException:
        return None
    except Exception:
        return None
    finally:
        signal.alarm(0)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_excel", help="Excel file (.xlsx)")
    parser.add_argument("output_sdf", help="Output SDF file")
    parser.add_argument("--confs", type=int, default=50, help="Conformers per molecule")
    parser.add_argument("--cores", type=int, default=32, help="Number of CPU cores to use")
    args = parser.parse_args()

    print(f"Reading {args.input_excel}...")
    
    df = pd.read_excel(args.input_excel).dropna(how='all')
    df = df.drop_duplicates(subset=df.columns[0], keep='first')
    
    tasks = []
    for i in range(len(df)):
        try:
            # Explicitly force Column 0 as ID and Column 1 as SMILES
            mol_id = str(df.iloc[i, 0]).strip()
            smiles = str(df.iloc[i, 1]).strip()
            
            if mol_id.lower() != 'nan' and smiles.lower() != 'nan':
                tasks.append((smiles, mol_id, args.confs))
        except IndexError:
            continue

    print(f"Generating conformers for {len(tasks)} inputs using {args.cores} cores...")
    
    with multiprocessing.Pool(args.cores) as pool:
        results = pool.map(process_molecule, tasks)

    print("Writing valid conformers to SDF with guaranteed IDs...")
    writer = Chem.SDWriter(args.output_sdf)
    valid_mols = 0
    total_poses = 0
    
    for result in results:
        if result is not None:
            mol, mol_id = result
            if mol.GetNumConformers() > 0:
                # STAMP THE ID ON THE MOLECULE IMMEDIATELY BEFORE WRITING
                # This guarantees the top line of the SDF block will contain your original ID
                mol.SetProp("_Name", mol_id)
                
                for conf in mol.GetConformers():
                    writer.write(mol, confId=conf.GetId())
                    total_poses += 1
                valid_mols += 1
            
    writer.close()
    
    print(f"Successfully processed {valid_mols} molecules, generating {total_poses} total 3D poses.")

if __name__ == "__main__":
    main()