import os
import json
import csv
import argparse

def extract_scores(boltz_out_dir):
    """Recursively finds confidence JSONs and extracts the ipTM score."""
    results = {}
    
    # Iterate over each main job folder (e.g., Rags_WT, Rags_W1A)
    for folder_name in os.listdir(boltz_out_dir):
        folder_path = os.path.join(boltz_out_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue
            
        # Walk recursively through all subdirectories inside this job folder
        found_json = False
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                if file.startswith("confidence_") and file.endswith(".json"):
                    filepath = os.path.join(root, file)
                    try:
                        with open(filepath, 'r') as f:
                            data = json.load(f)
                            # Extract the ipTM score
                            iptm = data.get('iptm', 0.0) 
                            results[folder_name] = iptm
                            found_json = True
                    except Exception as e:
                        print(f"Error reading {filepath}: {e}")
                    break # Stop looking at files in this specific directory
            
            if found_json:
                break # Stop walking subdirectories for this job once found
                
    return results

def main():
    parser = argparse.ArgumentParser(description="Parse nested Boltz-2 JSON outputs and calculate delta ipTM.")
    parser.add_argument("-d", "--dir", required=True, help="Directory containing Boltz-2 output folders.")
    parser.add_argument("-o", "--output", default="alanine_scan_results.csv", help="Output CSV filename.")
    parser.add_argument("-w", "--wildtype", required=True, help="A string to identify the Wild-Type folder (e.g., '_WT').")
    
    args = parser.parse_args()
    
    print(f"Scanning directory: {args.dir}...")
    scores = extract_scores(args.dir)
    
    if not scores:
        print("No confidence JSON files found. Check your directory path.")
        return

    # 1. Identify the Wild-Type score
    wt_key = None
    for key in scores.keys():
        if args.wildtype in key:
            wt_key = key
            break
            
    if not wt_key:
        print(f"Error: Could not find a folder containing '{args.wildtype}'.")
        print(f"Available folders: {list(scores.keys())}")
        return
        
    wt_score = scores[wt_key]
    print(f"Wild-Type ({wt_key}) ipTM: {wt_score:.4f}")
    
    # 2. Calculate the Delta ipTM for all mutants
    final_data = []
    for mutant, score in scores.items():
        if mutant == wt_key:
            continue
        
        delta_iptm = score - wt_score
        final_data.append({
            'Mutant': mutant,
            'ipTM': score,
            'Delta_ipTM': delta_iptm
        })
        
    # 3. Sort by the largest drop in score (most negative Delta_ipTM)
    final_data.sort(key=lambda x: x['Delta_ipTM'])
    
    # 4. Write to CSV
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Mutant', 'ipTM', 'Delta_ipTM'])
        writer.writeheader()
        
        writer.writerow({'Mutant': wt_key, 'ipTM': wt_score, 'Delta_ipTM': 0.0})
        writer.writerows(final_data)
        
    print(f"\nDone! Results sorted by highest impact and saved to {args.output}")

if __name__ == "__main__":
    main()