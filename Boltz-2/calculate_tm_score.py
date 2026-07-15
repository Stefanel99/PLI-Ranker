#!/usr/bin/env python3
from Bio.PDB import MMCIFParser, Superimposer
import numpy as np
from pathlib import Path
import json
from tqdm import tqdm

def get_ca_coords(cif_path, chain_id=None):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    model = next(structure.get_models())
    
    if chain_id is None:
        chain_id = list(model.child_dict.keys())[0]
    
    chain = model[chain_id]
    ca_atoms = [res["CA"] for res in chain if "CA" in res]
    
    if len(ca_atoms) == 0:
        raise ValueError(f"No CA atoms found in chain '{chain_id}' of {cif_path}")
    
    return ca_atoms, chain_id

def tm_score(native_cif, model_cif, native_chain=None, model_chain=None):
    native_ca, native_chain_id = get_ca_coords(native_cif, native_chain)
    model_ca, model_chain_id = get_ca_coords(model_cif, model_chain)
    
    L_native = len(native_ca)
    L_model = len(model_ca)
    L_aligned = min(L_native, L_model)
    
    native_ca = native_ca[:L_aligned]
    model_ca = model_ca[:L_aligned]
    
    if L_aligned < 10:
        raise ValueError(f"Too few aligned residues ({L_aligned}) for TM-score")
    
    sup = Superimposer()
    sup.set_atoms(native_ca, model_ca)
    sup.apply(model_ca)
    
    native_coords = np.array([a.get_coord() for a in native_ca])
    model_coords = np.array([a.get_coord() for a in model_ca])
    dists = np.linalg.norm(native_coords - model_coords, axis=1)
    
    L_target = len(native_ca)
    D0 = 1.24 * (L_target - 15) ** (1.0 / 3.0) - 1.8
    D0 = max(D0, 0.5)
    
    terms = 1.0 / (1.0 + (dists / D0) ** 2)
    tm_score_val = terms.sum() / float(L_target)
    
    pct_aligned = float(L_aligned) / float(L_native) * 100.0
    
    return float(tm_score_val), native_chain_id, model_chain_id, int(L_aligned), int(L_native), float(pct_aligned)

def find_complexes(reference_dir, predictions_dir):
    ref_files = sorted(list(Path(reference_dir).glob("*.cif")))
    results = []
    missing_preds = []
    
    print(f"\nReference files found: {len(ref_files)}")
    
    for ref_file in ref_files:
        complex_name = ref_file.stem
        pdb_code = complex_name.split("_")[0]
        pred_dir = Path(predictions_dir) / pdb_code
        
        print(f"{complex_name} -> looking in {pred_dir}")
        
        if not pred_dir.exists():
            missing_preds.append(complex_name)
            print(f"WARNING: Prediction directory not found: {pred_dir}")
            continue
        
        pred_files = sorted(list(pred_dir.glob(f"*{complex_name}*.cif")))
        
        if len(pred_files) == 0:
            missing_preds.append(complex_name)
            print(f"WARNING: No prediction files matching '{complex_name}'")
            continue
        
        print(f"Found {len(pred_files)} predictions")
        results.append({
            "complex": complex_name,
            "native": str(ref_file),
            "predictions": [str(p) for p in pred_files]
        })
    
    if missing_preds:
        print(f"\nMISSING predictions for {len(missing_preds)} complexes:")
        for m in missing_preds:
            print(f" - {m}")
    
    return results

def main():
    reference_dir = "reference_cifs"
    predictions_dir = "boltz2_results/cif_files"
    output_file = "tm_scores_results.json"
    
    complexes = find_complexes(reference_dir, predictions_dir)
    
    print(f"\nProcessing {len(complexes)} complexes with predictions")
    
    all_results = []
    
    for cx in tqdm(complexes, desc="Processing complexes"):
        native_file = cx["native"]
        
        for pred_file in cx["predictions"]:
            try:
                score, n_chain, m_chain, L_aligned, L_native, pct_aligned = tm_score(native_file, pred_file)
                
                all_results.append({
                    "complex": cx["complex"],
                    "prediction_file": Path(pred_file).name,
                    "tm_score": score,
                    "native_chain": n_chain,
                    "model_chain": m_chain,
                    "residues_aligned": L_aligned,
                    "residues_native": L_native,
                    "pct_aligned": round(pct_aligned, 2)
                })
            except Exception as e:
                all_results.append({
                    "complex": cx["complex"],
                    "prediction_file": Path(pred_file).name,
                    "tm_score": None,
                    "error": str(e)
                })
    
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\nSaved {len(all_results)} results to {output_file}")
    
    valid_scores = [r["tm_score"] for r in all_results if r["tm_score"] is not None]
    if valid_scores:
        print(f"Mean TM-score: {np.mean(valid_scores):.4f}")
        print(f"Median TM-score: {np.median(valid_scores):.4f}")
    
    valid_pct = [r["pct_aligned"] for r in all_results if r.get("pct_aligned") is not None]
    if valid_pct:
        print(f"Mean alignment coverage: {np.mean(valid_pct):.1f}%")

if __name__ == "__main__":
    main()
