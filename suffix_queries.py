import sys
from rdkit import Chem

def main():
    # Update this to your actual ZINC conformer SDF file
    input_sdf = "East_face2_top_hits.sdf" 
    output_sdf = "East_face2_top_hits_suffixed.sdf"
    
    print(f"Reading {input_sdf} and suffixing conformer names...")
    
    supplier = Chem.ForwardSDMolSupplier(input_sdf, removeHs=False)
    writer = Chem.SDWriter(output_sdf)
    
    name_counts = {}
    total_poses = 0
    
    for mol in supplier:
        if mol is None: 
            continue
            
        base_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "ZINC_Unknown"
        
        # Track how many times we've seen this base ID and append the count
        if base_name not in name_counts:
            name_counts[base_name] = 1
        else:
            name_counts[base_name] += 1
            
        unique_name = f"{base_name}_{name_counts[base_name]}"
        mol.SetProp("_Name", unique_name)
        
        writer.write(mol)
        total_poses += 1

    writer.close()
    
    print("\n--- Suffixing Complete ---")
    print(f"Total poses renamed: {total_poses}")
    print(f"Unique parent molecules: {len(name_counts)}")
    print(f"Saved to {output_sdf}. You are ready for alignment!")

if __name__ == "__main__":
    main()