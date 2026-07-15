#!/usr/bin/env python3
import json
import sys
import httpx
import asyncio
import logging
import os
import yaml
from fastapi import HTTPException
from typing import Any, Dict
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{task_id}"
PUBLIC_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

async def make_nvcf_call(function_url: str, data: Dict[str, Any]):
    api_key = os.getenv("NGC_API_KEY")
    if not api_key:
        raise ValueError("Error: NGC_API_KEY environment variable is not set.")
    
    max_retries = 5
    for attempt in range(max_retries):
        async with httpx.AsyncClient() as client:
            headers = {
                "Authorization": f"Bearer {api_key}",
                "NVCF-POLL-SECONDS": "300",
                "Content-Type": "application/json",
                "Accept": "application/json"
            }
            response = await client.post(function_url, json=data, headers=headers, timeout=400)
            
            if response.status_code == 429:
                wait_time = (2 ** attempt) + 5
                print(f"Rate limited. Retrying in {wait_time} seconds...")
                await asyncio.sleep(wait_time)
                continue
            if response.status_code == 202:
                task_id = response.headers.get("nvcf-reqid")
                while True:
                    await asyncio.sleep(5)
                    status_res = await client.get(STATUS_URL.format(task_id=task_id), headers=headers, timeout=400)
                    if status_res.status_code == 200:
                        logger.info("== Task Completed ==")
                        exec_details = status_res.headers.get("x-nvcf-execution-details")
                        if exec_details:
                            logger.info(f"Execution details header: {exec_details}")
                        else:
                            logger.info("No 'x-vcf-execution-details' header found...")
                        try:
                            full_json = status_res.json()
                            logger.info(f"Full status response keys: {full_json.keys()}")
                            if 'responseHeader' in full_json:
                                logger.info(f"Response header metadata: {full_json['responseHeader']}")
                        except Exception as e:
                            logger.error(f"Could not parse JSON: {e}")
                        return status_res.status_code, status_res
                    elif status_res.status_code >= 400:
                        raise HTTPException(status_res.status_code, status_res.text)
            elif response.status_code == 200:
                return response.status_code, response
            else:
                raise HTTPException(status_code=response.status_code, detail=response.text)
    raise Exception("Max retries exceeded for 429 rate limit!")

async def main():
    if len(sys.argv) < 2:
        print("Usage: python run_boltz2.py input.yaml")
        return
    yaml_input = sys.argv[1]
    with open(yaml_input, 'r') as file:
        input_data = yaml.safe_load(file)
        
    seeds = [321, 394, 422, 948]
    for complex_data in input_data['complexes']:
        protein_id, ligand_id = complex_data['pdb_id'], complex_data['ligand_id']
        protein_seq = complex_data['protein_sequence']
        ligand_smile = complex_data['ligand_smile']
        
        if protein_seq and ligand_smile:
            for seed in seeds:
                print(f"Running the process for seed number: {seed}")
                output_dir = Path("boltz2_results/eval_param/")
                output_dir.mkdir(parents=True, exist_ok=True)
                file_path = output_dir / f"{protein_id}_{ligand_id}_s{seed}_results.json"
                
                if file_path.exists():
                    print(f"Skipping {file_path.name}, already exists!")
                    continue
                    
                data = {
                    "polymers": [
                        {
                            "id": protein_id,
                            "molecule_type": "protein",
                            "sequence": protein_seq,
                        },
                    ],
                    "ligands": [
                        {
                            "smiles": ligand_smile,
                            "id": ligand_id,
                        },
                    ],
                    "recycling_steps": 3,
                    "sampling_steps": 75,
                    "diffusion_samples": 25,
                    "step_scale": 1.0
                }
                
                print(f"Sending {protein_id} with {ligand_id} request to NVIDIA Cloud...")
                try:
                    code, response = await make_nvcf_call(PUBLIC_URL, data)
                    if code == 200:
                        json_data = response.json()
                        file_path.write_text(json.dumps(json_data, indent=4))
                        print(f"Done! Saved to {file_path}")
                        
                        structures = json_data.get("structures", [])
                        cif_output_dir = Path(f"boltz2_results/cif_files/{protein_id}")
                        cif_output_dir.mkdir(parents=True, exist_ok=True)
                        
                        if structures:
                            for rank, structure_object in enumerate(structures):
                                cif_content = structure_object.get("structure")
                                if cif_content:
                                    sample_cif_path = cif_output_dir / f"{protein_id}_{ligand_id}_s{seed}_rank{rank}.cif"
                                    sample_cif_path.write_text(cif_content)
                                    print(f"Saved CIF file: {sample_cif_path.name}")
                                else:
                                    print(f"No cif content found in structure rank {rank}")
                            print(f"Done! Saved JSON and CIF file for {protein_id}_{ligand_id} (SEED {seed}, RANK {rank})")
                        else:
                            print(f"Done! Saved JSON, but no structures found for {protein_id}")
                except Exception as e:
                    print(f"Failed to process {protein_id} with seed {seed}: {e}")
        else:
            print(f"Skipping {protein_id} due to missing data!")
    print(f"The seeds used in this case were: {seeds}")

if __name__ == "__main__":
    asyncio.run(main())
