import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog

def main():
    parser = argparse.ArgumentParser(description="Filter PAINS compounds from a ROSHAMBO2 alignment scores CSV.")
    parser.add_argument("--scores", required=True, help="Input CSV file containing the alignment scores")
    parser.add_argument("--library", required=True, help="Input Excel file containing the SMILES library")
    parser.add_argument("--output", default="PAINS_free_scores.csv", help="Output filename for the clean results")
    args = parser.parse_args()

    print(f"Loading alignment scores from {args.scores}...")
    try:
        scores_df = pd.read_csv(args.scores)
    except FileNotFoundError:
        print(f"Error: Could not find scores file: {args.scores}")
        return

    print(f"Loading SMILES library from {args.library}...")
    try:
        lib_df = pd.read_excel(args.library).dropna(how='all')
    except FileNotFoundError:
        print(f"Error: Could not find library file: {args.library}")
        return
        
    lib_df = lib_df.drop_duplicates(subset=lib_df.columns[0], keep='first')
    
    # Create a fast lookup dictionary: { 'ID': 'SMILES' }
    id_to_smiles = dict(zip(
        lib_df.iloc[:, 0].astype(str).str.strip(), 
        lib_df.iloc[:, 1].astype(str).str.strip()
    ))

    print("Initializing RDKit PAINS filters...")
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)

    # Automatically detect the target column name
    target_col = 'Name' if 'Name' in scores_df.columns else 'Target_ID'
    if target_col not in scores_df.columns:
        target_col = scores_df.columns[1] # Fallback to the second column

    def is_pains_free(target_id):
        target_id = str(target_id).strip()
        smiles = id_to_smiles.get(target_id)
        
        if not smiles or smiles.lower() == 'nan':
            return False  
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False  
            
        # Return True only if it does NOT match a PAINS alert
        return not catalog.HasMatch(mol)

    print("Screening all candidates for structural alerts...")
    initial_count = len(scores_df)
    
    mask = scores_df[target_col].apply(is_pains_free)
    clean_df = scores_df[mask]
    
    dropped_count = initial_count - len(clean_df)

    clean_df.to_csv(args.output, index=False)
    
    print("\n--- Triage Complete ---")
    print(f"Total alignment pairs evaluated: {initial_count}")
    print(f"PAINS artifacts removed:         {dropped_count}")
    print(f"Clean, assay-ready candidates:   {len(clean_df)}")
    print(f"Saved highly confident list to:  {args.output}")

if __name__ == "__main__":
    main()