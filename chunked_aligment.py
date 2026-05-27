import os
import sys
import argparse
import pandas as pd
from rdkit import Chem
from roshambo2 import Roshambo2

def main():
    parser = argparse.ArgumentParser(description="Chunked GPU alignment using ROSHAMBO2 to prevent memory overflow.")
    parser.add_argument("--query", required=True, help="Input SDF file containing the 3D queries (e.g., zinc_queries_suffixed.sdf)")
    parser.add_argument("--target", required=True, help="Target HDF5 file containing the in-house library (e.g., ETP_cpds_clean.h5)")
    parser.add_argument("--output", default="master_scores.csv", help="Final stitched CSV output filename")
    parser.add_argument("--chunk_size", type=int, default=150, help="Number of queries to process per batch (default: 150)")
    
    args = parser.parse_args()

    print(f"Reading {args.query} and splitting into chunks of {args.chunk_size}...")
    
    try:
        supplier = Chem.ForwardSDMolSupplier(args.query, removeHs=False)
    except OSError:
        print(f"Error: Could not find or open input file '{args.query}'")
        sys.exit(1)
        
    chunk_idx = 1
    current_chunk_mols = []
    all_results = []

    for mol in supplier:
        if mol is not None:
            current_chunk_mols.append(mol)
            
        if len(current_chunk_mols) == args.chunk_size:
            process_chunk(current_chunk_mols, chunk_idx, args.target, all_results)
            chunk_idx += 1
            current_chunk_mols = []
            
    # Process any remaining poses at the end of the file
    if current_chunk_mols:
        process_chunk(current_chunk_mols, chunk_idx, args.target, all_results)

    # Prevent the concatenation crash if absolutely every chunk failed
    if not all_results:
        print("\nFatal Error: No chunks completed successfully. Nothing to stitch.")
        sys.exit(1)

    print(f"\nStitching {len(all_results)} chunks together...")
    
    # Because of our dictionary-wrapping fix in process_chunk, this will succeed
    final_df = pd.concat(all_results, ignore_index=True)
    
    # Save the master results file
    final_df.to_csv(args.output, index=False)
    print(f"SUCCESS! Master scores saved to {args.output}")

def process_chunk(mols, chunk_idx, target_h5, all_results):
    # Using a unique temp file name prevents race conditions
    chunk_filename = f"temp_query_chunk_{chunk_idx}.sdf"
    
    # Write the small batch to disk for ROSHAMBO to read
    writer = Chem.SDWriter(chunk_filename)
    for m in mols:
        writer.write(m)
    writer.close()
    
    print(f"\n--- Launching GPU Alignment for Chunk {chunk_idx} ({len(mols)} poses) ---")
    try:
        aligner = Roshambo2(chunk_filename, target_h5, color=True)
        
        # 1. Compute the raw dictionary on the GPU
        raw_results = aligner.compute(backend='cuda', n_gpus=1)
        
        # 2. Bulletproof pandas conversion for scalar dictionaries
        try:
            scores_df = pd.DataFrame(raw_results)
        except ValueError:
            # Wrap the dictionary in a list if pandas throws the scalar/index error
            scores_df = pd.DataFrame([raw_results])
        
        # 3. Append the clean DataFrame
        all_results.append(scores_df)
        
    except Exception as e:
        print(f"Error during chunk {chunk_idx}: {e}")
    finally:
        # Clean up the temp file
        if os.path.exists(chunk_filename):
            os.remove(chunk_filename)

if __name__ == "__main__":
    main()