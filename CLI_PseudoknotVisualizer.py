from coloring import CLI_coloring_canonical, load_colors_from_json, get_color_for_depth
from argparser import argparser, args_validation
from analysis.parsers import raw_df_processing, filter_abnormal_pairs
from config import RNAVIEW_DIR, RNAVIEW_EXEC, PseudoKnotVisualizer_DIR, INTERMEDIATE_DIR, DSSR_EXEC
from rna import PKextractor
from addressRNAviewOutput import load_rnaview_data #, extract_base_pairs_from_rnaview,
from addressDSSROutput import load_dssr_data #, extract_base_pairs_from_dssr,
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.MMCIFParser import MMCIFParser
import subprocess
import os
import pathlib
import shutil

colors = load_colors_from_json(PseudoKnotVisualizer_DIR / "colors.json")
# rnaview_exec = RNAVIEW_EXEC

class _ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_model(self, model):
        return True
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    def accept_residue(self, residue):
        return True
    def accept_atom(self, atom):
        return True


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

    # チェーン限定のPDBを書き出して RNAView に渡す（チェーン混入を避ける）
    arg = "--pdb"
    print(f"rnaview starts with {struct_file} and chain {chain_id}, output type is {arg} (chain-scoped PDB)")
    copied_file = pathlib.Path(INTERMEDIATE_DIR) / f"{pathlib.Path(struct_file).stem}_chain_{chain_id}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(copied_file), _ChainSelect(chain_id))

    subprocess.run(
        [RNAVIEW_EXEC, "-p", arg, str(copied_file)],
        env={"RNAVIEW": RNAVIEW_DIR},
        cwd=INTERMEDIATE_DIR,
        check=True
        )

    print("rnaview done.")
    result_file = pathlib.Path(INTERMEDIATE_DIR) / (copied_file.name + ".out")
    df = load_rnaview_data(str(result_file))
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
    if result.returncode != 0 or not json_output_path.exists():
        print(f"DSSR failed with return code: {result.returncode}")
        print(f"stdout: {result.stdout}")
        print(f"stderr: {result.stderr}")
    return load_dssr_data(str(json_output_path))

def CLI_PseudoKnotVisualizer(pdb_file, chain_id, format, output_file, model_id, annotator="RNAView", include_all=False):
    # パーサーの選択に応じてベースペアを抽出
    if annotator.upper() == "DSSR":
        raw_df = CLI_dssr(pdb_file, chain_id)
    
    elif annotator.upper() == "RNAVIEW":
        raw_df = CLI_rnaview(pdb_file, chain_id)
    # チェーンフィルタ（対象チェーン内のペアのみ残す）
    try:
        if not raw_df.empty and "chain1" in raw_df.columns and "chain2" in raw_df.columns:
            before = len(raw_df)
            filtered = raw_df[(raw_df["chain1"] == chain_id) & (raw_df["chain2"] == chain_id)].copy()
            if len(filtered) == 0 and before > 0:
                # フィルタで全て消えた場合は未フィルタのまま進める（チェーン限定PDBを使っているため）
                print(f"[CLI] Chain filter '{chain_id}' removed all pairs; keeping unfiltered {before} pairs")
            else:
                raw_df = filtered
                print(f"[CLI] Filtered by chain '{chain_id}': {len(raw_df)}/{before} pairs")
    except Exception as e:
        print(f"[CLI] Chain filtering skipped due to error: {e}")
    processed_df = raw_df_processing(raw_df, annotator)
    # remove abnormal pairs
    processed_df, abnormal_pairs, dup_canonical_pairs = filter_abnormal_pairs(processed_df)
    # Canonical-only by default unless include_all=True
    if not include_all:
        before_cnt = len(processed_df)
        processed_df = processed_df[processed_df["is_canonical"]].copy()
        print(f"[CLI] Using canonical base pairs only: {len(processed_df)}/{before_cnt}")
    else:
        canon_cnt = len(processed_df[processed_df["is_canonical"]])
        noncanon_cnt = len(processed_df) - canon_cnt
        print(f"[CLI] Using all base pairs: {len(processed_df)} total ({canon_cnt} canonical, {noncanon_cnt} non-canonical)")

    # print(f"Processed DataFrame:\n{processed_df.head()}")
    # Orientation normalization: ensure (i < j)
    BPL = []
    for _, row in processed_df.iterrows():
        i, j = row["position"]
        BPL.append((i, j) if i < j else (j, i))
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    PKlayers = PKextractor(BPL)

    with open(output_file, "w") as f:
        # 1) Precoloring (whiten target first)
        if format.lower() == "pymol":
            f.write(f"color white, {pdb_id} and chain {chain_id}\n")
        elif format.lower() == "chimera":
            # Chimera chain-wide whitening (model required)
            if model_id is None:
                print("[CLI] Warning: Chimera format requested but model_id is None; skipping whitening.")
            else:
                f.write(f"color white #{model_id}:.{chain_id}\n")

        # 2) Layer coloring + selection commands
        for depth, PKlayer in enumerate(PKlayers):
            color = get_color_for_depth(depth + 1, colors)
            script = CLI_coloring_canonical(pdb_id, model_id, chain_id, PKlayer, color, format)
            f.write(script)

            # Add paper-friendly selections for PyMOL output
            if format.lower() == "pymol":
                all_res = []
                for i, j in PKlayer:
                    all_res.extend([str(i), str(j)])
                if all_res:
                    res_expr = "+".join(all_res)
                    paper_name = "core" if depth == 0 else f"pk{depth}"
                    paper_name = f"{str(pathlib.Path(pdb_file).stem)}_{chain_id}_{paper_name}"
                    f.write(f"select {paper_name}, {pdb_id} and chain {chain_id} and resi {res_expr}\n")

    print("Coloring done.")
    print(f"Depth is {len(PKlayers)}")
    print(f"Output script is saved as {output_file}")
    return

def main():
    args = argparser()
    args_validation(args)

    print("PseudoKnotVisualizer started.")
    # Backward compatibility: accept legacy --parser if present
    annotator = getattr(args, 'annotator', None) or getattr(args, 'parser', 'RNAView')
    CLI_PseudoKnotVisualizer(args.input, args.chain, args.format, args.output, args.model, annotator, include_all=getattr(args, 'include_all', False))
    print("PseudoKnotVisualizer finished: " + args.output)

if __name__ == "__main__":
    main()
