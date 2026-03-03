import argparse
import os

def extract_sequences(filepath):
    """Extracts amino acid sequences directly from SEQRES or ATOM lines."""
    sequences = {}
    
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    with open(filepath, 'r') as f:
        for line in f:
            # 1. Try to read from SEQRES lines first (Fastest and cleanest)
            if line.startswith("SEQRES"):
                parts = line.split()
                if len(parts) >= 5:
                    chain_id = parts[2]
                    residues = parts[4:] 
                    
                    if chain_id not in sequences:
                        sequences[chain_id] = []
                        
                    for res in residues:
                        if res in aa_map:
                            sequences[chain_id].append(aa_map[res])

            # 2. Fallback: If no SEQRES lines exist, read from ATOM lines
            elif line.startswith("ATOM  ") and " CA " in line:
                chain_id = line[21].strip()
                resname = line[17:20].strip()
                
                if not chain_id:
                    chain_id = "UNKNOWN"
                    
                if chain_id not in sequences:
                    sequences[chain_id] = []
                    
                if resname in aa_map:
                    sequences[chain_id].append(aa_map[resname])

    final_sequences = {chain: "".join(seq) for chain, seq in sequences.items()}
    if not final_sequences:
        print("Warning: Could not find any valid sequences in the file.")
    return final_sequences

def parse_residue_ranges(range_str):
    """Parses a string like '[(1-20),25,27,35,(50-67)]' into a list of integers."""
    clean_str = range_str.replace(' ', '').replace('[', '').replace(']', '').replace('(', '').replace(')', '')
    indices = set()
    
    if not clean_str:
        return sorted(list(indices))
        
    parts = clean_str.split(',')
    for part in parts:
        if '-' in part:
            start, end = part.split('-')
            indices.update(range(int(start), int(end) + 1))
        else:
            indices.add(int(part))
            
    return sorted(list(indices))

def write_yaml(sequences, output_path):
    """Writes a Boltz-2 formatted YAML file using the updated version 1 schema."""
    lines = [
        "version: 1",
        "sequences:"
    ]
    for chain_id, seq in sequences.items():
        lines.append(f"  - protein:")
        lines.append(f"      id: {chain_id}")
        lines.append(f"      sequence: {seq}")
        
    with open(output_path, 'w') as f:
        f.write("\n".join(lines) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Prepare Boltz-2 YAML files for targeted Alanine Scanning.")
    parser.add_argument("-i", "--input", required=True, help="Input PDB or mmCIF file.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for YAML files.")
    parser.add_argument("-b", "--basename", required=True, help="Base name for the output YAML files.")
    parser.add_argument("-c", "--chain", required=True, help="Target chain ID to mutate (e.g., B).")
    parser.add_argument("-r", "--residues", required=False, 
                        help="Residues to mutate (1-based index). Format: '1-20,25,27,35,50-67'. If omitted, scans the whole chain.")
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    sequences = extract_sequences(args.input)
    
    if args.chain not in sequences:
        raise ValueError(f"Chain '{args.chain}' not found. Available chains: {list(sequences.keys())}")
        
    # 1. Write Wild-Type YAML
    wt_path = os.path.join(args.outdir, f"{args.basename}_WT.yaml")
    write_yaml(sequences, wt_path)
    print(f"Generated Wild-Type: {wt_path}")
    
    # 2. Determine target indices (1-based sequence indexing)
    target_seq = sequences[args.chain]
    if args.residues:
        target_indices = parse_residue_ranges(args.residues)
    else:
        target_indices = list(range(1, len(target_seq) + 1))
        
    # 3. Write Alanine Mutants
    for idx in target_indices:
        i = idx - 1 # Convert 1-based index to 0-based Python string index
        
        if i < 0 or i >= len(target_seq):
            print(f"Warning: Requested residue index {idx} is out of bounds for chain {args.chain} (length {len(target_seq)}). Skipping.")
            continue
            
        res = target_seq[i]
        if res == 'A':
            print(f"Skipping index {idx}: Residue is already Alanine.")
            continue
            
        mut_seq = target_seq[:i] + 'A' + target_seq[i+1:]
        mut_sequences = sequences.copy()
        mut_sequences[args.chain] = mut_seq
        
        mut_name = f"{res}{idx}A"
        mut_path = os.path.join(args.outdir, f"{args.basename}_{args.chain}_{mut_name}.yaml")
        
        write_yaml(mut_sequences, mut_path)
        print(f"Generated Mutant: {mut_path}")

if __name__ == "__main__":
    main()