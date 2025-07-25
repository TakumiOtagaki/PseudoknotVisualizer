from config import RNAVIEW_DIR, RNAVIEW_EXEC, PseudoKnotVisualizer_DIR, INTERMEDIATE_DIR, DSSR_EXEC
from coloring import coloring_canonical, load_colors_from_json
from argparser import argparser, args_validation
from analysis.parsers import raw_df_processing, filter_abnormal_pairs
from rna import PKextractor
import os
from pymol import cmd
import tempfile
import subprocess

from addressRNAviewOutput import load_rnaview_data
from addressDSSROutput import load_dssr_data
import pathlib

# DEBUG = True
DEBUG = False

colors = load_colors_from_json(PseudoKnotVisualizer_DIR / "colors.json")

def clear_intermediate_files(except_files=[]):
    # intermediate dir には他のゴミのファイルがあるので消しておく
    for f in os.listdir(INTERMEDIATE_DIR):
        # .gitkeep は残す
        if f == ".gitkeep":
            continue
        if f.endswith(".out") or f.endswith(".pdb") or f.endswith(".ps") or f.endswith(".xml") or f.endswith(".cif") or f.endswith(".pdb_new"):
            if f not in except_files:
                # os.remove(INTERMEDIATE_DIR + f)
                os.remove(pathlib.Path(INTERMEDIATE_DIR) / f)
    return

def is_pure_rna(pdb_object, chain=None):
    """
    指定した pdb_object( + chain) が
    ・標準的な RNA 塩基 (A, C, G, U, I)
    ・あるいは水分子やイオンなど (HOH, MG, NA, K, CL, etc.)
    以外を一切含まない場合に True を返す関数。
    
    - pdb_object: PyMOL上にロードされているオブジェクト名 (例: "1ABC")
    - chain: チェーンID を文字列で指定 (例: "A"). 
                Noneの場合は全チェーンを対象にチェックする。
    """
    # ▼ ここに標準的な RNA の残基名を定義 (実際にはより多くの修飾などがある場合は追加)
    rna_bases = {"A", "C", "G", "U", "I"}
    
    # ▼ ここに「無視してOK」な水・イオンの残基名を列挙
    #   (HOH/WAT/H2O は水、MG, NA, CL, K, CA, ZN ...etc.)
    ignore_resn = {"HOH", "WAT", "H2O", "MG", "NA", "K", "CL", "CA", "ZN"}

    # チェーンを指定したい場合は selection を絞る
    selection = pdb_object
    if chain:
        selection = f"{pdb_object} and chain {chain}"

    # PyMOLの model を取得
    model = cmd.get_model(selection)

    for atom in model.atom:
        # アトムの持つ残基名 (resn) をチェック
        resn_upper = atom.resn.upper()  # 大文字にしておく
        if resn_upper in ignore_resn:
            # 水やイオンなどはスキップ
            continue
        if resn_upper not in rna_bases:
            # RNA 塩基以外が見つかった場合は False
            return False
    
    # ループを抜けた = 問題ある残基が見つからなかった場合
    return True

def check_residues_start_from_one(pdb_object, chain):
    """
    指定したpdb_objectの chain のレジデュー番号が1から始まっているかをチェックする。
    """
    model = cmd.get_model(f"{pdb_object} and chain {chain}")
    resi_list = []
    for atom in model.atom:
        try:
            r = int(atom.resi)
            if r not in resi_list:
                resi_list.append(r)
        except ValueError:
            raise ValueError(f"Invalid residue number format in {pdb_object} chain {chain}: {atom.resi}")
            pass
    
    if not resi_list:
        return False
    
    return min(resi_list) == 1

def auto_renumber_residues(pdb_object, chain):
    """
    指定したpdb_objectの chain に対して、
    レジデュー番号の最小値が1でない場合に、自動で 1 から始まるように再番号付けする。
    例: min_resiが100の場合、 100 -> 1, 101 -> 2, ... のようにalterする。
    """
    # pdb_object and chain のアトム情報を取得
    model = cmd.get_model(f"{pdb_object} and chain {chain}")

    # ユニークな整数レジデュー番号のみを収集 (insertion codeなどは一旦無視)
    resi_list = []
    for atom in model.atom:
        try:
            r = int(atom.resi)
            if r not in resi_list:
                resi_list.append(r)
        except ValueError:
            # resi が変な文字付き(例: 101A 等)の場合、スキップ or 追加処理要検討
            pass

    if not resi_list:
        # 見つからなかった場合はスキップ
        print(f"[auto_renumber_residues] {pdb_object}, chain {chain} でレジデュー番号が見つかりませんでした。")
        return

    min_resi = min(resi_list)
    
    # すでに1から始まっていれば何もしない
    if min_resi == 1: return
    
    # (min_resi - 1) を引くことによって、最小値を1に合わせる
    offset = min_resi - 1
    
    print(f"[auto_renumber_residues] renumbering chain {chain} by offset={offset} (min_resi was {min_resi})")
    
    # alter で一括変更 (resiを int(resi)-offset に置き換える)
    cmd.alter(
        f"{pdb_object} and chain {chain}",
        f"resi=str(int(resi)-{offset})"
    )
    
    # チェーン間の干渉を避けるため、全体のソートは削除
    # cmd.sort(pdb_object)

def rnaview_wrapper(pdb_object, chain):
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTERMEDIATE_DIR) as tmp_pdb:
            pdb_path = tmp_pdb.name # tmp.pdb is created and deleted automatically after the block.
            cmd.save(pdb_path, f"{pdb_object} and chain {chain}", format="pdb")

            result = subprocess.run(
                [RNAVIEW_EXEC, "-p", "--pdb", pdb_path],
                env={"RNAVIEW": RNAVIEW_DIR},
                cwd=INTERMEDIATE_DIR,
                check=True
            )
            if result.returncode != 0:
                raise Exception("RNAVIEW failed")
    except Exception as e:
        raise Exception("RNAVIEW failed or Exporting PDB failed: " + str(e))
    result_file = pathlib.Path(INTERMEDIATE_DIR) / (pathlib.Path(pdb_path).name + ".out")
    raw_df = load_rnaview_data(result_file)
    return raw_df


def dssr_wrapper(pdb_object, chain):
    """DSSR wrapper function to extract base pairs"""
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTERMEDIATE_DIR) as tmp_pdb:
            pdb_path = tmp_pdb.name
            cmd.save(pdb_path, f"{pdb_object} and chain {chain}", format="pdb")

            # DSSR実行（JSONフォーマットで出力）
            json_output_path = pathlib.Path(INTERMEDIATE_DIR) / (pathlib.Path(pdb_path).name + ".dssr.json")
            result = subprocess.run(
                [str(DSSR_EXEC), f"-i={pdb_path}", "--json", f"-o={json_output_path}"],
                cwd=INTERMEDIATE_DIR,
                check=True,
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                raise Exception("DSSR failed")
    except Exception as e:
        raise Exception("DSSR failed or Exporting PDB failed: " + str(e))

    raw_df = load_dssr_data(json_output_path)
    return raw_df


def PseudoKnotVisualizer(pdb_object, chain=None, parser="RNAView", auto_renumber=True, only_pure_rna=False, skip_precoloring=False, selection=True):
    """
    PseudoKnotVisualizer: Visualizing Pseudo Knots in RNA structure.
    Usage: pkv pdb_object, chain=A, parser=RNAView
     - pdb_object(str): PDB object name
     - chain(str) : Chain ID of the RNA structure.
        If not specified, all chains will be analyzed.
     - parser(str) [default: "RNAView"]: Structure parser to use ("DSSR" or "RNAView").
     - auto_renumber(bool) [auto_renumber: True]: If True, automatically renumber residues from 1,
        to avoid the error caused by non-sequential residue numbers in the input PDB file.
     - only_pure_rna(bool) [default: False]: If True, only standard RNA bases (A, C, G, U, I) are analyzed.
     - skip_precoloring(bool) [default: False]: If True, all atoms are not colored 'white' before coloring the base pairs.
     - selection(bool): If True, selection will be created for each layer: pdb_object_pkorder0, pdb_object_pkorder1, pdb_object_pkorder2, ...
    """
    print("version ", PseudoKnotVisualizer_DIR / "VERSION.txt")
    
    if chain is None:
        chains = cmd.get_chains(pdb_object)
        print("Chain ID is not specified and there are multiple chains. All chains ID will be analyzed: " + ", ".join(chains))
        for chain in chains:
            PseudoKnotVisualizer(pdb_object, chain, auto_renumber, only_pure_rna, skip_precoloring, selection, parser)
        return
    elif chain not in cmd.get_chains(pdb_object):
        print(f"Chain {chain} is not found in the pdb object.")
        print(f"Available chains are: {', '.join(cmd.get_chains(pdb_object))}")
        return
    # ★ RNAViewを使用する場合のみ、レジデュー番号をチェックして必要に応じて補正
    if auto_renumber and parser.upper() == "RNAVIEW":
        if not check_residues_start_from_one(pdb_object, chain):
            print(f"[PseudoKnotVisualizer] Chain {chain}: residue numbers do not start from 1. Auto-renumbering will be applied.")
            auto_renumber_residues(pdb_object, chain)
        else:
            # print(f"[PseudoKnotVisualizer] Chain {chain}: residue numbers start from 1. No auto-renumbering needed.")
            pass
    if only_pure_rna:
        if not is_pure_rna(pdb_object, chain):
            print("The structure contains non-standard RNA bases or other molecules.")
            print("If you want to analyze them, please set pure_rna=False.")
            return
    
    # パーサーの選択に応じてベースペアを抽出
    if parser.upper() == "DSSR":
        raw_df = dssr_wrapper(pdb_object, chain)
    elif parser.upper() == "RNAVIEW":
        raw_df = rnaview_wrapper(pdb_object, chain)
    else:
        raise ValueError(f"Unsupported parser: {parser}. Use 'DSSR' or 'RNAView'.")
    
    processed_df = raw_df_processing(raw_df, parser)
    processed_df, abnormal_pairs, dup_canonical_pairs = filter_abnormal_pairs(processed_df)
    # print(processed_df)
    BPL = [tuple(row["position"]) for _, row in processed_df.iterrows()]
    # print(f"extracted base pairs: {BPL}")
    PKlayers = PKextractor(BPL)

    if not skip_precoloring:
        print(f"Precoloring all atoms to white since skip_precoloring is {skip_precoloring}")
        cmd.color("white", f"{pdb_object} and chain {chain}")
    for depth, PKlayer in enumerate(PKlayers):
        color = colors[str(depth + 1)]
        print(f"Coloring layer {depth + 1} with color: {color}")
        
        # PKlayerから全ての残基番号を取得して、PyMOL用の選択文字列を作成
        all_residues = []
        for i, j in PKlayer:
            all_residues.extend([str(i), str(j)])
        selection_str = "+".join(all_residues)
        
        for i, j in PKlayer:
            coloring_canonical(pdb_object, chain, i, color)
            coloring_canonical(pdb_object, chain, j, color)
        if selection:
            # selectコマンドで選択オブジェクトを作成
            selection_name = f"{pdb_object}_l{depth}"
            cmd.select(selection_name, f"{pdb_object} and chain {chain} and resi {selection_str}")
            print(f"Created selection: {selection_name} with residues {selection_str}")
        print(f"Layer {depth + 1}: (i, j) = {PKlayer}")
    print("Coloring done.")
    print(f"Depth is {len(PKlayers)}")
    
    clear_intermediate_files()
    return


cmd.extend("PseudoKnotVisualizer", PseudoKnotVisualizer)
cmd.extend("pkv", PseudoKnotVisualizer)

if __name__ == "__main__":
    pass

