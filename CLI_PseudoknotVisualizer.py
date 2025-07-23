from PseudoknotVisualizer import clear_intermediate_files, colors
from coloring import CLI_coloring_canonical, load_colors_from_json
from argparser import argparser, args_validation
from analysis.parsers import raw_df_processing, filter_abnormal_pairs
from config import RNAVIEW_DIR, RNAVIEW_EXEC, PseudoKnotVisualizer_DIR, INTERMEDIATE_DIR, DSSR_EXEC
from rna import PKextractor
from addressRNAviewOutput import load_rnaview_data #, extract_base_pairs_from_rnaview,
from addressDSSROutput import load_dssr_data #, extract_base_pairs_from_dssr,
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
import tempfile
import subprocess
import os
import pathlib
import shutil

colors = load_colors_from_json(PseudoKnotVisualizer_DIR / "colors.json")
# rnaview_exec = RNAVIEW_EXEC

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
    copied_file = pathlib.Path(INTERMEDIATE_DIR) / pathlib.Path(struct_file).name
    shutil.copy2(struct_file, copied_file)


    result = subprocess.run(
        [RNAVIEW_EXEC, "-p", arg, copied_file],
        env={"RNAVIEW": RNAVIEW_DIR},
        cwd=INTERMEDIATE_DIR,
        check=True
        )
    if result.returncode != 0:
        raise Exception("RNAVIEW failed")

    print("rnaview done.")
    result_file = pathlib.Path(INTERMEDIATE_DIR) / (pathlib.Path(struct_file).name + ".out")
    df = load_rnaview_data(result_file)
    # valid_bps_df = extract_canonicalbp_from_rnaview(df)
    # print(valid_bps_df)
    # BPL = [(row["left_idx"], row["right_idx"]) for _, row in valid_bps_df.iterrows()]
    # BPL = [(row["left_idx"], row["right_idx"]) for _, row in df.iterrows()]

    # return BPL
    return df

def CLI_dssr(struct_file, chain_id):
    """CLI version of DSSR wrapper"""
    # intermediate 以下に複製する
    copied_file = pathlib.Path(INTERMEDIATE_DIR) / pathlib.Path(struct_file).name
    shutil.copy2(struct_file, copied_file)

    print(f"DSSR starts with {struct_file}")
    
    # DSSR実行（JSONフォーマットで出力）
    json_output_path = pathlib.Path(INTERMEDIATE_DIR) / (pathlib.Path(struct_file).name + ".dssr.json")
    result = subprocess.run(
        [str(DSSR_EXEC), f"-i={str(copied_file)}", "--json", f"-o={str(json_output_path)}"],
        cwd=INTERMEDIATE_DIR,
        check=True,
        capture_output=True,
        text=True
    )
    print(f"command: \n {str(DSSR_EXEC)} -i={str(copied_file)} --json -o={str(json_output_path)}")
    # if result.returncode != 0:
    #     raise Exception("DSSR failed")
        # エラー内容を出力
    if result.returncode != 0 or not json_output_path.exists():
        print(f"DSSR failed with return code: {result.returncode}")
        print(f"stdout: {result.stdout}")
        print(f"stderr: {result.stderr}")
        return load_dssr_data(json_output_path)
        # raise Exception("DSSR failed")
    # print(f"DSSR output generated: {json_output_path}")

    print("DSSR done.")
    df = load_dssr_data(json_output_path)
    # valid_bps_df = extract_canonicalbp_from_dssr(json_output_path)
    # print(valid_bps_df)
    # BPL = [(row["left_idx"], row["right_idx"]) for _, row in valid_bps_df.iterrows()]

    return df

def CLI_PseudoKnotVisualizer(pdb_file, chain_id, format, output_file, model_id, parser="RNAView"):
    # パーサーの選択に応じてベースペアを抽出
    if parser.upper() == "DSSR":
        raw_df = CLI_dssr(pdb_file, chain_id)
    
    elif parser.upper() == "RNAVIEW":
        raw_df = CLI_rnaview(pdb_file, chain_id)
    else:
        raise ValueError(f"Unsupported parser: {parser}. Use 'DSSR' or 'RNAView'.")

    processed_df = raw_df_processing(raw_df, parser)
    # remove abnormal pairs
    processed_df, abnormal_pairs = filter_abnormal_pairs(processed_df)
    BPL = [(row["left_idx"], row["right_idx"]) for _, row in processed_df.iterrows()]
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
    
    clear_intermediate_files()
    return

def main():
    args = argparser()
    args_validation(args)

    print("PseudoKnotVisualizer started.")
    CLI_PseudoKnotVisualizer(args.input, args.chain, args.format, args.output, args.model, args.parser)
    print("PseudoKnotVisualizer finished: " + args.output)

if __name__ == "__main__":
    main()