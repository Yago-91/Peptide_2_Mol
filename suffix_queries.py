import argparse
import sys
from rdkit import Chem

def main():
    parser = argparse.ArgumentParser(description="Suffix conformer names in an SDF file to guarantee unique IDs for ROSHAMBO2.")
    parser.add_argument("--input", required=True, help="Input SDF file containing the 3D queries")
    parser.add_argument("--output", default="queries_suffixed.sdf", help="Output SDF file to save the renamed queries")
    args = parser.parse_args()

    print(f"Reading {args.input} and suffixing conformer names...")
    
    try:
        supplier = Chem.ForwardSDMolSupplier(args.input, removeHs=False)
    except OSError:
        print(f"Error: Could not find or open input file '{args.input}'")
        sys.exit(1)
        
    try:
        writer = Chem.SDWriter(args.output)
    except OSError:
        print(f"Error: Could not create output file '{args.output}'")
        sys.exit(1)
    
    name_counts = {}
    total_poses = 0
    
    for mol in supplier:
        if mol is None: 
            continue
            
        # Extract the base ID or assign a placeholder if missing
        base_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
        
        # Track how many times we've seen this base ID and append the incremental count
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
    print(f"Total poses renamed:       {total_poses}")
    print(f"Unique parent molecules:   {len(name_counts)}")
    print(f"Saved highly-unique file:  {args.output}")

if __name__ == "__main__":
    main()