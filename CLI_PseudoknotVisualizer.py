from PseudoknotVisualizer import clear_intermediate_files, rnaview, colors
from coloring import CLI_coloring_canonical, load_colors_from_json
from argparser import argparser, args_validation
from config import RNAVIEW_PATH, RNAVIEW, PseudoKnotVisualizer_DIR, INTEREMEDIATE_DIR
from rna import PKextractor
from addressRNAviewOutput import extract_base_pairs_from_rnaview
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
import tempfile
import subprocess
import os
import pathlib
import shutil

rnaview = RNAVIEW_PATH / "rnaview"

def CLI_rnaview(struct_file, chain_id):
    # 入力ファイルが .cif か .pdb かを拡張子で判定
    ext = os.path.splitext(struct_file)[1].lower()
    file_type = None
    if ext == ".cif" or ext == ".mmcif":
        parser = MMCIFParser(QUIET=True)
        file_type = "cif"
    elif ext == ".pdb":
        parser = PDBParser(QUIET=True)
        file_type = "pdb"
    else:
        raise ValueError("Input file should be .cif or .pdb")

    structure = parser.get_structure("structure", struct_file)
    chain_structure = structure[0]
    selected_chain = None

    for chain in chain_structure:
        if chain.id == chain_id:
            selected_chain = chain
            break

    if selected_chain is None:
        raise ValueError(f"Chain ID {chain_id} not found in {struct_file}")

    # PDB と CIF で出力方法を分ける
    arg = "--cif" if file_type == "cif" else "--pdb"

    print(f"rnaview starts with {struct_file} and chain {chain_id}, output type is {arg}")
    # intermediate 以下に複製する
    copied_file = pathlib.Path(INTEREMEDIATE_DIR) / pathlib.Path(struct_file).name
    shutil.copy2(struct_file, copied_file)


    result = subprocess.run(
        [rnaview, "-p", arg, copied_file],
        env={"RNAVIEW": RNAVIEW},
        cwd=INTEREMEDIATE_DIR,
        check=True
        )
    if result.returncode != 0:
        raise Exception("RNAVIEW failed")

    print("rnaview done.")
    result_file = pathlib.Path(INTEREMEDIATE_DIR) / (pathlib.Path(struct_file).name + ".out")
    valid_bps_df = extract_base_pairs_from_rnaview(result_file)
    print(valid_bps_df)
    BPL = [(row["left_idx"], row["right_idx"]) for _, row in valid_bps_df.iterrows()]

    return BPL

def CLI_PseudoKnotVisualizer(pdb_file, chain_id, format, output_file, model_id):
    BPL = CLI_rnaview(pdb_file, chain_id)
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    PKlayers = PKextractor(BPL)

    with open(output_file, "w") as f:
        for depth, PKlayer in enumerate(PKlayers):
            color = colors[str(depth + 1)]
            script = CLI_coloring_canonical(pdb_id, model_id, chain_id, PKlayer, color, format)
            f.write(script)

    print("Coloring done.")
    print(f"Depth is {len(PKlayers)}")
    print(f"Output script is saved as {output_file}")
    
    # clear_intermediate_files()
    return

def main():
    args = argparser()
    args_validation(args)

    print("PseudoKnotVisualizer started.")
    CLI_PseudoKnotVisualizer(args.input, args.chain, args.format, args.output, args.model)

if __name__ == "__main__":
    main()