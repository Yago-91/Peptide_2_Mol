import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Rank drug pairs using thermodynamic ensemble scoring.")
    parser.add_argument("--input", required=True, help="Input PAINS-free CSV")
    parser.add_argument("--output", default="Ensemble_Ranked_Hits.csv", help="Output ranked CSV")
    parser.add_argument("--threshold", type=float, default=1.0, help="Baseline score to be considered a 'good' pose")
    args = parser.parse_args()

    print(f"Loading data from {args.input}...")
    df = pd.read_csv(args.input)

    # 1. Strip the pose suffix (e.g., "ZINC123_42" -> "ZINC123")
    df['ZINC_Parent'] = df['Query_ID'].astype(str).str.rsplit('_', n=1).str[0]
    
    # 2. Dynamic column mapping for ROSHAMBO2 V2 outputs
    target_col = 'name' if 'name' in df.columns else df.columns[1]
    
    if 'tanimoto_combo_legacy' in df.columns:
        score_col = 'tanimoto_combo_legacy'
    else:
        print("Fatal Error: Could not find the combination score column.")
        return

    print(f"Using '{score_col}' for scoring and '{target_col}' for targets.")
    print("Calculating Ensemble metrics per drug pair...")

    # Define a custom aggregation function to calculate our math
    def calculate_ensemble_metrics(group):
        scores = group[score_col].sort_values(ascending=False).values
        
        max_score = scores[0] if len(scores) > 0 else 0
        
        # Metric 1: Count of poses above threshold
        good_poses = scores[scores > args.threshold]
        hit_count = len(good_poses)
        
        # Metric 2: Mean of Top 3 poses (Smooths out 1-hit wonders)
        top_3_mean = np.mean(scores[:3]) if len(scores) >= 3 else np.mean(scores)
        
        # Metric 3: Integrated Ensemble Score (Area over threshold)
        ies = np.sum(good_poses - args.threshold) if hit_count > 0 else 0
        
        return pd.Series({
            f'Max_{score_col}': max_score,
            'Top_3_Mean': top_3_mean,
            f'Poses_Over_{args.threshold}': hit_count,
            'Integrated_Ensemble_Score': ies,
            'Total_Poses_Evaluated': len(scores)
        })

    # Group by the unique pair (ZINC Parent + In-House Compound)
    grouped = df.groupby(['ZINC_Parent', target_col]).apply(calculate_ensemble_metrics).reset_index()

    # Sort by the Integrated Ensemble Score to reward high frequencies of high scores
    ranked_df = grouped.sort_values(by='Integrated_Ensemble_Score', ascending=False)

    ranked_df.to_csv(args.output, index=False)
    
    print("\n--- Ranking Complete ---")
    print(f"Total unique drug pairs evaluated: {len(ranked_df)}")
    print(f"Saved ensemble rankings to: {args.output}")

if __name__ == "__main__":
    main()