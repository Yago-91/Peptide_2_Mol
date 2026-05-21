import sys
from rdkit import Chem

def main():
    input_sdf = "ETP_cpds.sdf"
    output_sdf = "ETP_cpds_clean.sdf"
    
    print(f"Scanning {input_sdf} for degenerate 1-atom artifacts...")
    
    # Use a ForwardSDMolSupplier to stream the 9.2GB file without loading it all into RAM
    supplier = Chem.ForwardSDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    
    total_poses = 0
    written_poses = 0
    dropped_mols = set()

    for mol in supplier:
        if mol is None:
            continue
        
        total_poses += 1
        mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
        
        # A valid drug molecule must have more than 3 atoms to calculate a 3D covariance matrix
        if mol.GetNumAtoms() <= 3:
            dropped_mols.add(mol_id)
            continue
            
        writer.write(mol)
        written_poses += 1
        
        if total_poses % 500000 == 0:
            print(f"Processed {total_poses} poses...")

    writer.close()
    
    print("\n--- Cleaning Complete ---")
    print(f"Total poses scanned: {total_poses}")
    print(f"Valid poses retained: {written_poses}")
    print(f"Degenerate poses dropped: {total_poses - written_poses}")
    if dropped_mols:
        print(f"Dropped molecule IDs (fewer than 4 atoms): {', '.join(list(dropped_mols)[:10])}...")

if __name__ == "__main__":
    main()