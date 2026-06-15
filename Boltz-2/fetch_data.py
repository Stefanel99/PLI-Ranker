#!/usr/bin/env python3

import httpx
import re
import sys
import yaml
import json
from pathlib import Path


def download_ground_truth(protein_id:str,ligand_id:str,output_dir:str="reference_cifs"):
    output_path=Path(output_dir)
    output_path.mkdir(parents=True,exist_ok=True)
    save_path=output_path/f"{protein_id}_{ligand_id}.cif"
    if save_path.exists():
        print(f"File {save_path} already exists. Skipping download.")
        return str(save_path)