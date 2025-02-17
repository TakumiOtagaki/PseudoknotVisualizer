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

from pymol import cmd

from pymol import cmd

def is_pure_rna(pdb_object, chain_id=None):
    """
    指定した pdb_object(＋chain_id) が
    ・標準的な RNA 塩基 (A, C, G, U, I)
    ・あるいは水分子やイオンなど (HOH, MG, NA, K, CL, etc.)
    以外を一切含まない場合に True を返す関数。
    
    - pdb_object: PyMOL上にロードされているオブジェクト名 (例: "1ABC")
    - chain_id: チェーンID を文字列で指定 (例: "A"). 
                Noneの場合は全チェーンを対象にチェックする。
    """
    # ▼ ここに標準的な RNA の残基名を定義 (実際にはより多くの修飾などがある場合は追加)
    rna_bases = {"A", "C", "G", "U", "I"}
    
    # ▼ ここに「無視してOK」な水・イオンの残基名を列挙
    #   (HOH/WAT/H2O は水、MG, NA, CL, K, CA, ZN ...etc.)
    ignore_resn = {"HOH", "WAT", "H2O", "MG", "NA", "K", "CL", "CA", "ZN"}

    # チェーンを指定したい場合は selection を絞る
    selection = pdb_object
    if chain_id:
        selection = f"{pdb_object} and chain {chain_id}"

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

def auto_renumber_residues(pdb_object, chain_id):
    """
    指定したpdb_objectの chain_id に対して、
    レジデュー番号の最小値が1でない場合に、自動で 1 から始まるように再番号付けする。
    例: min_resiが100の場合、 100 -> 1, 101 -> 2, ... のようにalterする。
    """
    # pdb_object and chain_id のアトム情報を取得
    model = cmd.get_model(f"{pdb_object} and chain {chain_id}")

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
        print(f"[auto_renumber_residues] {pdb_object}, chain {chain_id} でレジデュー番号が見つかりませんでした。")
        return

    min_resi = min(resi_list)
    
    # すでに1から始まっていれば何もしない
    if min_resi == 1: return
    
    # (min_resi - 1) を引くことによって、最小値を1に合わせる
    offset = min_resi - 1
    
    print(f"[auto_renumber_residues] renumbering chain {chain_id} by offset={offset} (min_resi was {min_resi})")
    
    # alter で一括変更 (resiを int(resi)-offset に置き換える)
    cmd.alter(
        f"{pdb_object} and chain {chain_id}",
        f"resi=str(int(resi)-{offset})"
    )
    
    # sort して整合性を取っておく (resi 順序などが正しく並ぶように)
    cmd.sort(pdb_object)

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


def PseudoKnotVisualizer(pdb_object, chain_id=None, auto_renumber=True, only_pure_rna=False):
    """
    PseudoKnotVisualizer: Visualizing Pseudo Knots in RNA structure.
    Usage: pkv pdb_object [,chain_id]
     - pdb_object(str): PDB object name
     - chain_id(str): Chain ID of the RNA structure.
        If not specified, all chains will be analyzed.
     - auto_renumber(bool): If True, automatically renumber residues from 1,
        to avoid the error caused by non-sequential residue numbers in the input PDB file.
     - only_pure_rna(bool): If True, only standard RNA bases (A, C, G, U, I) are analyzed.
    """
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
    # ★ ここで自動でレジデュー番号を補正する
    if auto_renumber:
        auto_renumber_residues(pdb_object, chain_id)
    if only_pure_rna:
        if not is_pure_rna(pdb_object, chain_id):
            print("The structure contains non-standard RNA bases or other molecules.")
            print("If you want to analyze them, please set pure_rna=False.")
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

