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

def CLI_rnaview(pdb_file, chain_id):
    # 入力ファイルが .cif か .pdb かを拡張子で判定
    ext = os.path.splitext(pdb_file)[1].lower()
    if ext == ".cif" or ext == ".mmcif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("structure", pdb_file)
    chain_structure = structure[0]
    selected_chain = None

    for chain in chain_structure:
        if chain.id == chain_id:
            selected_chain = chain
            break

    if selected_chain is None:
        raise ValueError(f"Chain ID {chain_id} not found in {pdb_file}")

    # PDB と CIF で出力方法を分ける
    if ext == ".cif" or ext == ".mmcif":
        io = MMCIFIO()
        suffix = ".cif"
        rnaview_args = [rnaview, "-p", "--cif"]  # CIF 用
    else:
        io = PDBIO()
        suffix = ".pdb"
        rnaview_args = [rnaview, "-p", "--pdb"]  # PDB 用

    io.set_structure(selected_chain)

    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix, dir=INTEREMEDIATE_DIR) as tmp_file:
            tmp_path = tmp_file.name
            io.save(tmp_path)

            # RNAVIEW 実行コマンドに一時ファイルを渡す
            rnaview_args.append(tmp_path)

            result = subprocess.run(
                rnaview_args,
                env={"RNAVIEW": RNAVIEW},
                cwd=INTEREMEDIATE_DIR,
                check=True
            )
            if result.returncode != 0:
                raise Exception("RNAVIEW failed")
    except Exception as e:
        raise Exception("RNAVIEW failed or Exporting structure failed: " + str(e))

    result_file = INTEREMEDIATE_DIR + os.path.basename(tmp_path) + ".out"
    valid_bps_df = extract_base_pairs(result_file)
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