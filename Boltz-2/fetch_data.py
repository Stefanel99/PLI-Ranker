#!/usr/bin/env python3

from urllib import response
from xmlrpc import client

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
    url=f"https://files.rcsb.org/view/{protein_id}.cif"
    try:
        with httpx.Client(follow_directions=True) as client:
            response=client.get(url,timeout=20)
            if response.status_code==200:
                save_path.write_text(response.text)
                print(f"Succesfully downloaded: {protein_id} !")
                return str(save_path)
            else:
                print(f"Failed to download {protein_id}. \nStatus code: {response.status_code}")
                return None
    except httpx.RequestError as e:
        print(f"Network error downloading {protein_id}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
    return None


def extract_smile_complex(ligand_id:str)->str:
    ligand_id=ligand_id.strip().upper()
    ligand_url=f"https://files.rcsb.org/ligands/view/{ligand_id}.cif"
    try:
        with httpx.Client(follow_redirects=True) as client:
            response=client.get(ligand_url,timeout=15)
        if response.status_code==200:
            clif_text=response.text
            for line in clif_text.splitlines():
                if "SMILES" in line and "CACTVS" in line:
                    match=re.search(r'"([^"]+)"',line)
                    if match :return match.group(1)
    except Exception as e:
        print(f" Error accessing RSCB for : {ligand_id}!")
        print(f" Failed to fetch SMILE for : {ligand_id} !")
    return False


def extract_protein_seq(protein_id:str)->str:
    protein_id=protein_id.strip().lower()
    protein_url=f"https://www.rcsb.org/fasta/entry/{protein_id}/display"
    try:
        with httpx.Client(follow_redirects=True) as client:
            response=client.get(protein_url,timeout=15)
            if response.status_code!=200:
                print(f"Error: Could not find PDB {protein_id} (Status: {response})")
                return None
            fasta_text=response.text.strip()
            if not fasta_text:
                print(f"Error: Received empty FASTA for {protein_id}")
                return None
            lines=fasta_text.splitlines()
            if len(lines)<2:
                print(f"Error: Malformed FASTA for {protein_id}")
                return None
            protein_seq="".join(lines.strip() for line in lines[1:])
            if not protein_seq:
                print(f"Error: No sequence data found in FASTA for {protein_id}")
                return None
            return protein_seq
    except httpx.ConnectTimeout as e:
        print(f"Error: Connection to RCSB timed out for {protein_id}!")
    except httpx.RequestError as e:
        print(f"Error: Network request failed for {protein_id}: {e}")
    except Exception as e:
        print(f"Unexpected error extracting sequence for {protein_id}!")
    return None



def fetch_data(input_json:str,output_yaml:str="complexes.yaml"):
    yaml_data={"complexes":[]}
    with open(input_json,"r") as file:
        data=json.load(file).get('complexes')
        for item in data:
            protein_id=item['pdb_id']
            ligand_id=item['ligand_id']
            download_ground_truth(protein_id,ligand_id)
            protein_seq=extract_protein_seq(protein_id)
            ligand_smile=extract_smile_complex(ligand_id)
            yaml_data['complexes'].append({
                "pdb_id":protein_id,
                "ligand_id":ligand_id,
                "protein_sequence":protein_seq,
                "ligand_smile":ligand_smile
            })
        output_path=Path(output_yaml)
        with output_path.open("w") as file:
            yaml.dump(yaml_data,file,sort_keys=False)
        print(f"Saved YAML to {output_path.resolve()}")



if __name__=="__main__":
    if len(sys.argv)<2:
        print("Usage: python fetch_data.py <input_json> [output_yaml]")
        sys.exit(1)
    input_file=sys.argv[1]
    output_file=sys.argv[2] if len(sys.argv)>2 else "complexes.yaml"
    fetch_data(input_file,output_file)