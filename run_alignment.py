import argparse
from roshambo2 import Roshambo2

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True, help="SDF file containing ZINC hits")
    parser.add_argument("--target", required=True, help="SDF file containing in-house library")
    parser.add_argument("--output_prefix", default="roshambo_results", help="Prefix for output files")
    args = parser.parse_args()

    print(f"Loading queries from {args.query} and targets from {args.target}...")
    
    # Initialize the ROSHAMBO2 engine
    # backend='cuda' pushes the volumetric integrals to the GPU
    # color=True calculates pharmacophoric feature overlap alongside steric shape
    aligner = Roshambo2(
        query_file=args.query, 
        target_file=args.target,
        backend='cuda',
        color=True
    )

    print("Executing GPU-accelerated alignment...")
    # This computes the ComboTanimoto (Shape + Color) for all pairs
    scores_df = aligner.compute()
    
    # Save the scores to a CSV for filtering and analysis
    csv_out = f"{args.output_prefix}_scores.csv"
    scores_df.to_csv(csv_out, index=False)
    print(f"Scores saved to {csv_out}")

    print("Extracting best fit aligned structures...")
    # Writes the physically aligned 3D coordinates of the best matches to a new SDF
    # You can open this file in PyMOL alongside your ZINC hits to visually inspect the overlap
    aligner.write_best_fit_structures(hits_sdf_prefix=args.output_prefix)
    print("Alignment complete. Best fits saved as SDF.")

if __name__ == "__main__":
    main()