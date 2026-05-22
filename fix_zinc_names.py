from rdkit import Chem

def main():
    input_sdf = "Filtered_hits_smiles.sdf"
    output_sdf = "zinc_unique.sdf"
    
    print(f"Scanning {input_sdf} for duplicate IDs...")
    
    supplier = Chem.ForwardSDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    
    seen_names = set()
    current_name = None
    actual_name_to_write = None
    suffix_counter = 2
    
    duplicates_fixed = 0

    for mol in supplier:
        if mol is None: 
            continue
            
        raw_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
        
        # Detect if we have moved to a new block of conformers
        if raw_name != current_name:
            current_name = raw_name
            
            # If we have already seen this ID earlier in the file, it's a duplicate row
            if current_name in seen_names:
                while f"{current_name}_v{suffix_counter}" in seen_names:
                    suffix_counter += 1
                actual_name_to_write = f"{current_name}_v{suffix_counter}"
                duplicates_fixed += 1
            else:
                actual_name_to_write = current_name
                suffix_counter = 2
                
            seen_names.add(actual_name_to_write)
            
        # Rename the molecule to the unique version
        mol.SetProp("_Name", actual_name_to_write)
        writer.write(mol)

    writer.close()
    print(f"Complete! Renamed {duplicates_fixed} duplicate groups.")
    print(f"Saved unique queries to {output_sdf}")

if __name__ == "__main__":
    main()