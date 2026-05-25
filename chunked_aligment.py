import os
import pandas as pd
from rdkit import Chem
from roshambo2 import Roshambo2

def main():
    input_query_sdf = "zinc_queries_suffixed.sdf"
    target_h5 = "ETP_cpds_clean.h5"
    final_output = "ZINC_ETP_master_scores.csv"
    
    # 150 queries against 1M targets limits matrix size to ~12GB RAM/VRAM
    chunk_size = 150 

    print(f"Reading {input_query_sdf} and splitting into GPU-safe chunks...")
    supplier = Chem.ForwardSDMolSupplier(input_query_sdf, removeHs=False)
    
    chunk_idx = 1
    current_chunk_mols = []
    all_results = []

    for mol in supplier:
        if mol is not None:
            current_chunk_mols.append(mol)
            
        if len(current_chunk_mols) == chunk_size:
            process_chunk(current_chunk_mols, chunk_idx, target_h5, all_results)
            chunk_idx += 1
            current_chunk_mols = []
            
    # Process the remaining poses at the end of the file
    if current_chunk_mols:
        process_chunk(current_chunk_mols, chunk_idx, target_h5, all_results)

    print(f"\nStitching {len(all_results)} chunks together...")
    final_df = pd.concat(all_results, ignore_index=True)
    
    # Save the massive results file
    final_df.to_csv(final_output, index=False)
    print(f"SUCCESS! Matrix limit bypassed. Master scores saved to {final_output}")

def process_chunk(mols, chunk_idx, target_h5, all_results):
    chunk_filename = f"temp_query_chunk.sdf"
    
    # Write the small batch to disk for ROSHAMBO to read
    writer = Chem.SDWriter(chunk_filename)
    for m in mols:
        writer.write(m)
    writer.close()
    
    print(f"\n--- Launching GPU Alignment for Chunk {chunk_idx} ({len(mols)} poses) ---")
    try:
        # Initialize engine on just this small chunk
        aligner = Roshambo2(chunk_filename, target_h5, color=True)
        scores_df = aligner.compute(backend='cuda', n_gpus=1)
        all_results.append(scores_df)
    except Exception as e:
        print(f"Error during chunk {chunk_idx}: {e}")
    finally:
        # Clean up the temp file
        if os.path.exists(chunk_filename):
            os.remove(chunk_filename)

if __name__ == "__main__":
    main()