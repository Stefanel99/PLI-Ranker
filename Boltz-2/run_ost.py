#!/usr/bin/env python3
import os
import subprocess
import glob
import shutil
import sys

if len(sys.argv) > 1:
    boltz_dir = sys.argv[1]
    ref_dir = sys.argv[2]
    output_dir = sys.argv[3]
else:
    ref_dir = "reference_cifs"
    boltz_dir = "boltz2_results/cif_files"
    output_dir = "ost_results"

os.makedirs(output_dir, exist_ok=True)

ost_path = shutil.which("ost")
if not ost_path:
    print("ERROR: 'ost' not found. Activate the conda env where ost is installed.")
    sys.exit(1)

print(f"Using ost: {ost_path}")

ref_files = glob.glob(os.path.join(ref_dir, "*.cif"))
print(f"Found {len(ref_files)} reference files\n")

for ref_path in ref_files:
    filename = os.path.basename(ref_path)
    pdb_id = filename.split("_")[0]
    target_subdir = os.path.join(boltz_dir, pdb_id)
    
    print(f"Processing {pdb_id}...")
    
    if not os.path.isdir(target_subdir):
        print("SKIP: Directory not found")
        continue
    
    cif_files = glob.glob(os.path.join(target_subdir, "*.cif"))
    if not cif_files:
        print("SKIP: No cif files found")
        continue
    
    print(f"Found {len(cif_files)} prediction files")
    
    for pred_path in cif_files:
        pred_filename = os.path.basename(pred_path)
        base_name = pred_filename.replace(".cif", "")
        output_file = os.path.join(output_dir, f"eval_{base_name}.json")
        
        cmd = [ost_path, "compare-ligand-structures", "-m", pred_path, "-r", ref_path, "--rmsd", "--fault-tolerant", "-o", output_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"OK: {pred_filename}")
        else:
            print(f"FAIL: {pred_filename} - {result.stderr.strip()}")

print("\nDone")
