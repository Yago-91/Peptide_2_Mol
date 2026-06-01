import argparse
import sys
import os
import pandas as pd
from multiprocessing import Pool

# Global variable inside each worker process to hold the catalog
worker_catalog = None

def init_worker():
    """Initializes the RDKit catalog once per CPU core when the process starts."""
    global worker_catalog
    from rdkit.Chem import FilterCatalog
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    worker_catalog = FilterCatalog.FilterCatalog(params)

def check_smiles_worker(smiles_val):
    """The parallel worker function executing on an individual CPU core."""
    from rdkit import Chem
    smiles = str(smiles_val).strip()
    
    if not smiles or smiles.lower() == 'nan' or smiles == '':
        return False
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
        
    # Return True if it does NOT match a PAINS alert
    return not worker_catalog.HasMatch(mol)

def main():
    parser = argparse.ArgumentParser(description="Parallel PAINS filtering utilizing multi-core CPUs.")
    parser.add_argument("--input", required=True, help="Input master scores CSV file")
    parser.add_argument("--output", default="PAINS_free_scores.csv", help="Output filename")
    parser.add_argument("--cores", type=int, default=None, help="Number of CPU cores to use (defaults to all available)")
    args = parser.parse_args()

    print(f"Loading alignment scores from {args.input}...")
    try:
        df = pd.read_csv(args.input)
    except FileNotFoundError:
        print(f"Error: Could not find input file: {args.input}")
        sys.exit(1)

    total_rows = len(df)
    if total_rows == 0:
        print("Error: Input file is empty.")
        sys.exit(1)

    if len(df.columns) < 3:
        print(f"Error: Input file lacks a 3rd column for SMILES. Found {len(df.columns)} columns.")
        sys.exit(1)

    smiles_col_name = df.columns[2]
    
    # Determine CPU core allocation
    available_cores = os.cpu_count()
    selected_cores = args.cores if args.cores else available_cores
    print(f"Detected SMILES in column: '{smiles_col_name}'")
    print(f"Allocating {selected_cores}/{available_cores} CPU cores for parallel screening...")

    mask = []
    processed_count = 0

    # Launch the parallel pool
    # The initializer ensures the heavy PAINS catalog is loaded only ONCE per core, not per SMILES
    with Pool(processes=selected_cores, initializer=init_worker) as pool:
        
        # imap maintains row order while processing asynchronously
        for result in pool.imap(check_smiles_worker, df[smiles_col_name], chunksize=100):
            mask.append(result)
            processed_count += 1
            
            # Live terminal progress updates
            if processed_count % 50 == 0 or processed_count == total_rows:
                percent = (processed_count / total_rows) * 100
                bar_length = 30
                filled_length = int(round(bar_length * processed_count / float(total_rows)))
                bar = '█' * filled_length + '-' * (bar_length - filled_length)
                
                sys.stdout.write(f"\rProgress: |{bar}| {percent:.1f}% ({processed_count}/{total_rows} rows)")
                sys.stdout.flush()

    print("\n\nFiltering dataset and writing output...")
    clean_df = df[mask]
    dropped_count = total_rows - len(clean_df)

    clean_df.to_csv(args.output, index=False)
    
    print("--- Triage Complete ---")
    print(f"Total alignment pairs evaluated: {total_rows}")
    print(f"PAINS artifacts removed:         {dropped_count}")
    print(f"Clean, assay-ready candidates:   {len(clean_df)}")
    print(f"Saved highly confident list to:  {args.output}")

if __name__ == "__main__":
    main()