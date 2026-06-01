import argparse
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog

def main():
    parser = argparse.ArgumentParser(description="Filter PAINS compounds using SMILES embedded directly in the scores CSV.")
    parser.add_argument("--input", required=True, help="Input master scores CSV file (e.g., ZINC_ETP_master_scores.csv)")
    parser.add_argument("--output", default="PAINS_free_scores.csv", help="Output filename for the clean results")
    args = parser.parse_args()

    print(f"Loading alignment scores from {args.input}...")
    try:
        df = pd.read_csv(args.input)
    except FileNotFoundError:
        print(f"Error: Could not find input file: {args.input}")
        sys.exit(1)

    if len(df.columns) < 3:
        print(f"Error: The input file does not have at least 3 columns. Found {len(df.columns)} columns.")
        sys.exit(1)

    # Identify the SMILES column dynamically (the 3rd column, which is index 2)
    smiles_col_name = df.columns[2]
    print(f"Detected SMILES data in the third column: '{smiles_col_name}'")

    print("Initializing RDKit PAINS filters...")
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)

    def is_pains_free(smiles_val):
        smiles = str(smiles_val).strip()
        
        if not smiles or smiles.lower() == 'nan' or smiles == '':
            return False  # Filter out rows with missing structural data
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False  # Filter out rows where RDKit cannot parse the chemistry
            
        # Return True only if it does NOT match a PAINS alert
        return not catalog.HasMatch(mol)

    print("Screening all rows for structural alerts...")
    initial_count = len(df)
    
    # Apply the filter directly using the third column
    mask = df[smiles_col_name].apply(is_pains_free)
    clean_df = df[mask]
    
    dropped_count = initial_count - len(clean_df)

    clean_df.to_csv(args.output, index=False)
    
    print("\n--- Triage Complete ---")
    print(f"Total alignment pairs evaluated: {initial_count}")
    print(f"PAINS artifacts removed:         {dropped_count}")
    print(f"Clean, assay-ready candidates:   {len(clean_df)}")
    print(f"Saved highly confident list to:  {args.output}")

if __name__ == "__main__":
    main()