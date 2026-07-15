#!/usr/bin/env python3
import json
import sys
import pandas as pd
import re
from pathlib import Path

def parse_boltz_id(file_path: Path):
    return file_path.stem.replace("_results", "")

def parse_ost_parts(file_path: Path):
    stem = file_path.stem
    parts = stem.split('_')
    if len(parts) >= 5:
        complex_id = f"{parts[1]}_{parts[2]}_{parts[3]}"
        rank_id = parts[4].replace("rank", "")
        return complex_id, rank_id
    return None, None

def extract_boltz_data(boltz_root: Path):
    data = []
    json_files = list(boltz_root.rglob("*.json"))
    for file_path in json_files:
        try:
            with open(file_path, 'r') as f:
                content = json.load(f)
            complex_id = parse_boltz_id(file_path)
            keys = [
                'confidence_scores', 'ptm_scores', 'iptm_scores',
                'ligand_iptm_scores', 'complex_plddt_scores',
                'complex_iplddt_scores', 'complex_pde_scores', 'complex_ipde_scores'
            ]
            available_arrays = [content.get(k) for k in keys if content.get(k) is not None]
            if not available_arrays: 
                continue
            max_len = len(available_arrays[0])
            for rank_idx in range(max_len):
                row = {'complex_id': str(complex_id), 'rank': str(rank_idx)}
                for k in keys:
                    arr = content.get(k, [])
                    row[f'boltz_{k}'] = arr[rank_idx] if rank_idx < len(arr) else None
                data.append(row)
        except Exception: 
            continue
    return pd.DataFrame(data)

def extract_ost_data(ost_root: Path):
    data = []
    json_files = list(ost_root.rglob("*.json"))
    for file_path in json_files:
        try:
            with open(file_path, 'r') as f:
                content = json.load(f)
            c_id, r_id = parse_ost_parts(file_path)
            if not c_id: 
                continue
            row = {'complex_id': str(c_id), 'rank': str(r_id)}
            assigned = content.get("rmsd", {}).get("assigned_scores", [])
            if assigned:
                row['ost_rmsd'] = assigned[0].get("bb_rmsd")
                row['ost_l_rmsd'] = assigned[0].get("score")
                row['ost_lddt_lp'] = assigned[0].get("lddt_lp")
                data.append(row)
        except Exception: 
            continue
    return pd.DataFrame(data)

def extract_tm_data(tm_file: Path):
    if not tm_file.exists():
        return pd.DataFrame()
    try:
        with open(tm_file, 'r') as f:
            content = json.load(f)
        data = []
        for entry in content:
            pred_file = entry.get('prediction_file', '')
            
            id_match = re.search(r'^(.*)\_rank', Path(pred_file).stem)
            full_complex_id = id_match.group(1) if id_match else entry.get('complex')
            
            rank_match = re.search(r'rank(\d+)', pred_file)
            rank_id = rank_match.group(1) if rank_match else "0"
            
            data.append({
                'complex_id': str(full_complex_id),
                'rank': str(rank_id),
                'tm_score': entry.get('tm_score'),
                'tm_coverage': entry.get('pct_aligned')
            })
        return pd.DataFrame(data)
    except Exception as e:
        print(f"Error parsing TM JSON: {e}")
        return pd.DataFrame()

def run_merge(b_path, o_path, tm_path, out_file):
    df_b = extract_boltz_data(Path(b_path))
    df_o = extract_ost_data(Path(o_path))
    df_tm = extract_tm_data(Path(tm_path))
    
    if df_b.empty or df_o.empty:
        print("Error: Boltz or OST data empty.")
        return
    
    merged = pd.merge(df_b, df_o, on=['complex_id', 'rank'], how='inner')
    
    if not df_tm.empty:
        df_tm['complex_id'] = df_tm['complex_id'].astype(str).str.strip()
        df_tm['rank'] = df_tm['rank'].astype(str).str.strip()
        merged['complex_id'] = merged['complex_id'].astype(str).str.strip()
        merged['rank'] = merged['rank'].astype(str).str.strip()
        
        merged = pd.merge(merged, df_tm, on=['complex_id', 'rank'], how='left')
    
    merged.to_csv(out_file, index=False)
    print(f"SUCCESS: {len(merged)} rows saved to {out_file}")

if __name__ == "__main__":
    run_merge(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
