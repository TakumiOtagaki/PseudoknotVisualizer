from config import RNAVIEW_PATH, RNAVIEW, WORK_DIR, INTEREMEDIATE_DIR
import os
from pymol import cmd
import tempfile
import subprocess

from addressRNAviewOutput import extract_base_pairs
# os.environ["RNAVIEW"] = RNAVIEW
# os.environ["RNAVIEW_PATH"] = RNAVIEW_PATH

def clear_intermediate_files():
    # intermediate dir には他のゴミのファイルがあるので消しておく
    for f in os.listdir(INTEREMEDIATE_DIR):
        if f.endswith(".out") or f.endswith(".pdb") or f.endswith(".ps") or f.endswith(".xml"):
            os.remove(INTEREMEDIATE_DIR + f)
    return

def rnaview(pdb_object, chain_id):
    chains = cmd.get_chains(pdb_object)
    if chain_id not in chains:
        raise Exception("Chain ID not found: should be one of " + ", ".join(chains))
    
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_pdb:
            pdb_path = tmp_pdb.name # tmp.pdb is created and deleted automatically after the block.
            cmd.save(pdb_path, pdb_object)

            result = subprocess.run(
                [RNAVIEW_PATH, pdb_path],
                env={"RNAVIEW": RNAVIEW},
                cwd=INTEREMEDIATE_DIR,
            )
            if result.returncode != 0:
                raise Exception("RNAVIEW failed")
    except Exception as e:
        raise Exception("RNAVIEW failed or Exporting PDB failed: " + str(e))

    result_file = INTEREMEDIATE_DIR + pdb_path.split("/")[-1] + ".out"
    valid_bps = extract_base_pairs(result_file) # pandas




        
        
        


    




