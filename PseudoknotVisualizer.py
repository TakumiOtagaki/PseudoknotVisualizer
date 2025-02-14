from config import RNAVIEW_PATH, RNAVIEW, PseudoKnotVisualizer_DIR, INTEREMEDIATE_DIR
from coloring import coloring_canonical, load_colors_from_json
from argparser import argparser, args_validation
from rna import PKextractor
import os
from pymol import cmd
import tempfile
import subprocess

from addressRNAviewOutput import extract_base_pairs
# os.environ["RNAVIEW"] = RNAVIEW
# os.environ["RNAVIEW_PATH"] = RNAVIEW_PATH

# DEBUG = True
DEBUG = False

colors = load_colors_from_json(PseudoKnotVisualizer_DIR + "/colors.json")

def clear_intermediate_files(except_files=[]):
    # intermediate dir には他のゴミのファイルがあるので消しておく
    for f in os.listdir(INTEREMEDIATE_DIR):
        if f.endswith(".out") or f.endswith(".pdb") or f.endswith(".ps") or f.endswith(".xml"):
            if f not in except_files:
                os.remove(INTEREMEDIATE_DIR + f)
    return

rnaview = os.path.join(RNAVIEW_PATH, "rnaview")
def rnaview_(pdb_object, chain_id):
    if DEBUG:
        pdb_id = pdb_object # for debugging
        chain_id = "A"
        pdb_path = "/large/otgk/PseudoknotVisualizer/intermediate/1KPD_test.pdb"
        result = subprocess.run(
            [rnaview, pdb_path],
            env={"RNAVIEW": RNAVIEW},
            cwd=INTEREMEDIATE_DIR,
            check=True
        )

    else:
        chains = cmd.get_chains(pdb_object)
        if chain_id not in chains:
            raise Exception("Chain ID not found: should be one of " + ", ".join(chains))
    
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTEREMEDIATE_DIR) as tmp_pdb:
                pdb_path = tmp_pdb.name # tmp.pdb is created and deleted automatically after the block.
                cmd.save(pdb_path, pdb_object)

                result = subprocess.run(
                    [rnaview, pdb_path],
                    env={"RNAVIEW": RNAVIEW},
                    cwd=INTEREMEDIATE_DIR,
                    check=True
                )
                if result.returncode != 0:
                    raise Exception("RNAVIEW failed")
        except Exception as e:
            raise Exception("RNAVIEW failed or Exporting PDB failed: " + str(e))

    result_file = INTEREMEDIATE_DIR + pdb_path.split("/")[-1] + ".out"
    valid_bps_df = extract_base_pairs(result_file) # pandas
    print(valid_bps_df)
    BPL = [(row["left_idx"], row["right_idx"]) for _, row in valid_bps_df.iterrows()]

    return BPL


def PseudoKnotVisualizer(pdb_object, chain_id):
    BPL = rnaview_(pdb_object, chain_id)
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
    # BPL = rnaview("1KPD", "A")
    # print(BPL)
    # PKlayers = PKextractor(BPL)
    # print(PKlayers)
    

    pass





        
        
        


    




