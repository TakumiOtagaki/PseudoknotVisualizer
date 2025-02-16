from config import RNAVIEW_PATH, RNAVIEW, PseudoKnotVisualizer_DIR, INTEREMEDIATE_DIR
from coloring import coloring_canonical, load_colors_from_json
from argparser import argparser, args_validation
from rna import PKextractor
import os
from pymol import cmd
import tempfile
import subprocess

from addressRNAviewOutput import extract_base_pairs_from_rnaview
import pathlib

# DEBUG = True
DEBUG = False

colors = load_colors_from_json(PseudoKnotVisualizer_DIR / "colors.json")
# rnaview = os.path.join(RNAVIEW_PATH, "rnaview")
rnaview = RNAVIEW_PATH / "rnaview"

def clear_intermediate_files(except_files=[]):
    # intermediate dir には他のゴミのファイルがあるので消しておく
    for f in os.listdir(INTEREMEDIATE_DIR):
        if f.endswith(".out") or f.endswith(".pdb") or f.endswith(".ps") or f.endswith(".xml") or f.endswith(".cif") or f.endswith(".pdb_new"):
            if f not in except_files:
                # os.remove(INTEREMEDIATE_DIR + f)
                os.remove(pathlib.Path(INTEREMEDIATE_DIR) / f)
    return

def rnaview_wrapper(pdb_object, chain_id):
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTEREMEDIATE_DIR) as tmp_pdb:
            pdb_path = tmp_pdb.name # tmp.pdb is created and deleted automatically after the block.
            cmd.save(pdb_path, pdb_object, format="pdb")

            result = subprocess.run(
                [rnaview, "-p", "--pdb", pdb_path],
                env={"RNAVIEW": RNAVIEW},
                cwd=INTEREMEDIATE_DIR,
                check=True
            )
            if result.returncode != 0:
                raise Exception("RNAVIEW failed")
    except Exception as e:
        raise Exception("RNAVIEW failed or Exporting PDB failed: " + str(e))

    # result_file = INTEREMEDIATE_DIR + pdb_path.split("/")[-1] + ".out"
    result_file = pathlib.Path(INTEREMEDIATE_DIR) / (pathlib.Path(pdb_path).name + ".out")
    valid_bps_df = extract_base_pairs_from_rnaview(result_file) # pandas
    print(valid_bps_df)
    BPL = [(row["left_idx"], row["right_idx"]) for _, row in valid_bps_df.iterrows()]

    return BPL


def PseudoKnotVisualizer(pdb_object, chain_id=None):
    if chain_id is None:
        chains = cmd.get_chains(pdb_object)
        print("Chain ID is not specified and there are multiple chains. All chains ID will be analyzed: " + ", ".join(chains))
        for chain_id in chains:
            PseudoKnotVisualizer(pdb_object, chain_id)
        return
    elif chain_id not in cmd.get_chains(pdb_object):
        print(f"Chain {chain_id} is not found in the pdb object.")
        print(f"Available chains are: {', '.join(cmd.get_chains(pdb_object))}")
        return
    BPL = rnaview_wrapper(pdb_object, chain_id)
    print(f"extracted base pairs: {BPL}")
    PKlayers = PKextractor(BPL)
    for depth, PKlayer in enumerate(PKlayers):
        # color = str(depth + 1)
        color = colors[str(depth + 1)]
        
        for i, j in PKlayer:
            coloring_canonical(pdb_object, chain_id, i, color)
            coloring_canonical(pdb_object, chain_id, j, color)
        print(f"Layer {depth + 1}: (i, j) = {PKlayer}")
    print("Coloring done.")
    print(f"Depth is {len(PKlayers)}")
    
    clear_intermediate_files()
    return


cmd.extend("PseudoKnotVisualizer", PseudoKnotVisualizer)
cmd.extend("pkv", PseudoKnotVisualizer)

if __name__ == "__main__":
    # clear_intermediate_files()
    # BPL = rnaview_wrapper("1KPD", "A")
    # print(BPL)
    # PKlayers = PKextractor(BPL)
    # print(PKlayers)
    pass

