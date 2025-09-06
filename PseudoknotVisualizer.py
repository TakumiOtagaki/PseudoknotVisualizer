from config import RNAVIEW_DIR, RNAVIEW_EXEC, PseudoKnotVisualizer_DIR, INTERMEDIATE_DIR, DSSR_EXEC
from coloring import coloring_canonical, load_colors_from_json, get_color_for_depth
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

def clear_intermediate_files(except_files=None):
    # intermediate dir には他のゴミのファイルがあるので消しておく
    if except_files is None:
        except_files = []
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
            # 挿入コード付きなどの非整数 resi は無視して判定する
            continue
    
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
        # print(f"[auto_renumber_residues] {pdb_object}, chain {chain} でレジデュー番号が見つかりませんでした。")
        print(f"[auto_renumber_residues] {pdb_object}, chain {chain} has not been found.")
        return

    min_resi = min(resi_list)
    
    # すでに1から始まっていれば何もしない
    if min_resi == 1: return
    
    # (min_resi - 1) を引くことによって、最小値を1に合わせる
    offset = min_resi - 1
    
    print(f"[auto_renumber_residues] renumbering chain {chain} by offset={offset} (min_resi was {min_resi})")
    
    # alter で一括変更 (resv を用いて安全に再番号付け)
    cmd.alter(
        f"{pdb_object} and chain {chain}",
        f"resi=str(resv-{offset})"
    )
    
    # チェーン間の干渉を避けるため、全体のソートは削除
    # cmd.sort(pdb_object)

def rnaview_wrapper(pdb_object, chain):
    try:
        # Check RNAView binary existence and guide user
        if not pathlib.Path(RNAVIEW_EXEC).exists():
            raise FileNotFoundError(
                "RNAView binary not found. Expected at:\n"
                f"  - {RNAVIEW_EXEC}\n\n"
                "How to fix:\n"
                "  1) Place the built 'rnaview' binary under 'RNAView/bin/rnaview' in this repository, or\n"
                "  2) Edit 'config.py' and set RNAVIEW_EXEC to your installation (e.g., '/opt/RNAView/bin/rnaview').\n"
                f"Also verify RNAVIEW_DIR points to: {RNAVIEW_DIR}"
            )
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTERMEDIATE_DIR) as tmp_pdb:
            pdb_path = tmp_pdb.name # tmp.pdb is created and deleted automatically after the block.
            cmd.save(pdb_path, f"{pdb_object} and chain {chain}", format="pdb")

            subprocess.run(
                [RNAVIEW_EXEC, "-p", "--pdb", pdb_path],
                env={"RNAVIEW": RNAVIEW_DIR},
                cwd=INTERMEDIATE_DIR,
                check=True
            )
    except Exception as e:
        raise Exception("RNAVIEW failed or Exporting PDB failed: " + str(e))
    result_file = pathlib.Path(INTERMEDIATE_DIR) / (pathlib.Path(pdb_path).name + ".out")
    raw_df = load_rnaview_data(str(result_file))
    return raw_df


def dssr_wrapper(pdb_object, chain):
    """DSSR wrapper function to extract base pairs"""
    try:
        # Check DSSR binary existence and guide user
        if not pathlib.Path(DSSR_EXEC).exists():
            raise FileNotFoundError(
                "DSSR binary not found. Expected at:\n"
                f"  - {DSSR_EXEC}\n\n"
                "How to fix:\n"
                "  1) Download 'x3dna-dssr' and place it under 'DSSR/x3dna-dssr' in this repository, or\n"
                "  2) Edit 'config.py' and set DSSR_EXEC to your installation (e.g., '/usr/local/bin/x3dna-dssr').\n"
                "On macOS, you may need: 'chmod +x x3dna-dssr' and allow it in System Settings > Privacy & Security."
            )
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=INTERMEDIATE_DIR) as tmp_pdb:
            pdb_path = tmp_pdb.name
            cmd.save(pdb_path, f"{pdb_object} and chain {chain}", format="pdb")

            # DSSR実行（JSONフォーマットで出力）
            json_output_path = pathlib.Path(INTERMEDIATE_DIR) / (pathlib.Path(pdb_path).name + ".dssr.json")
            subprocess.run(
                [str(DSSR_EXEC), f"-i={pdb_path}", "--json", f"-o={json_output_path}"],
                cwd=INTERMEDIATE_DIR,
                check=True,
                capture_output=True,
                text=True
            )
    except Exception as e:
        raise Exception("DSSR failed or Exporting PDB failed: " + str(e))

    raw_df = load_dssr_data(str(json_output_path))
    return raw_df


def PseudoKnotVisualizer(
    pdb_object,
    annotator="RNAView",
    chain=None,
    auto_renumber=True,
    only_pure_rna=False,
    skip_precoloring=False,
    selection=True,
    parser=None,  # deprecated: backward-compatible alias for annotator
    include_all=False,
):
    """
    PseudoKnotVisualizer: Visualize pseudoknot layers in RNA structures.

    PyMOL command:
        pkv object [,chain] [,annotator] [,auto_renumber] [,only_pure_rna] [,skip_precoloring] [,selection] [,include_all]

    Parameters
    ----------
    object : str
        Structure object name loaded in PyMOL.
    chain : str | None
        Chain ID. If omitted, all chains in the object are analyzed sequentially.
    annotator : {"RNAView", "DSSR"}
        Base-pair annotator. Default: "RNAView".
    auto_renumber : bool
        If True, renumber residues to start from 1 when necessary (mainly for RNAView) to avoid numbering issues.
    only_pure_rna : bool
        Reserved flag for future use. Currently NOT enforced and ignored at runtime.
    skip_precoloring : bool
        If True, do not pre-color the chain white before layer coloring.
    selection : bool
        If True, create selections per layer: "<obj>_c<chain>_l<depth>".
    include_all : bool
        If False (default), use canonical base pairs only (Watson-Crick + wobble). If True, include all pairs.

    Notes
    -----
    - Colors per layer are defined in 'colors.json'. If a layer index isn't configured,
      the 'default' color will be applied.
    - RNAView may require residues to start at 1. When annotator="RNAView" and auto_renumber=True,
      residues are renumbered per chain if necessary. For complex numbering, consider annotator="DSSR".
    - When include_all=True, non-canonical base pairs are included in addition to canonical ones.

    Deprecated
    ----------
    parser : str | None
        Old name for 'annotator'. Still accepted; if provided, it overrides 'annotator'.
    """

    # Backward compatibility: allow 'parser' keyword
    if parser is not None:
        print("[deprecated] 'parser' is deprecated. Use 'annotator' (\"RNAView\" or \"DSSR\").")
        annotator = parser
    # PyMOL からの引数が文字列の場合にも正しく解釈できるよう、bool へ正規化
    try:
        if not isinstance(include_all, bool):
            include_all = str(include_all).strip().lower() in ("1", "true", "t", "yes", "y", "on")
    except Exception:
        include_all = bool(include_all)
    # バージョン文字列の表示（存在しない場合は無視）
    try:
        with open(PseudoKnotVisualizer_DIR / "VERSION.txt", "r") as vf:
            version_str = vf.read().strip()
        print(f"version {version_str}")
    except Exception:
        pass
    print(
        f"arguments: pdb_object={pdb_object}, chain={chain}, annotator={annotator}, "
        f"auto_renumber={auto_renumber}, only_pure_rna={only_pure_rna}, skip_precoloring={skip_precoloring}, "
        f"selection={selection}, include_all={include_all}"
    )
    
    if chain is None:
        chains = cmd.get_chains(pdb_object)
        print("Chain ID is not specified and there are multiple chains. All chains ID will be analyzed: " + ", ".join(chains))
        for chain in chains:
            # Recurse per chain
            PseudoKnotVisualizer(pdb_object=pdb_object,
                                    chain=chain,
                                    auto_renumber=auto_renumber,
                                    only_pure_rna=only_pure_rna,
                                    skip_precoloring=skip_precoloring,
                                    selection=selection,
                                    annotator=annotator,
                                    include_all=include_all)
        return
    elif chain not in cmd.get_chains(pdb_object):
        print(f"Chain {chain} is not found in the pdb object.")
        print(f"Available chains are: {', '.join(cmd.get_chains(pdb_object))}")
        return
    # ★ RNAViewを使用する場合のみ、レジデュー番号をチェックして必要に応じて補正
    if auto_renumber and annotator.upper() == "RNAVIEW":
        if not check_residues_start_from_one(pdb_object, chain):
            # print(f"[PseudoKnotVisualizer] Chain {chain}: レジデュー番号が1から始まっていないため、RNAView用に補正します。")
            print(f"[PseudoKnotVisualizer] Chain {chain}: Residue numbers do not start from 1, renumbering for RNAView.")
            print("It may be better to use DSSR instead of RNAView.")
            auto_renumber_residues(pdb_object, chain)
        else:
            print(f"[PseudoKnotVisualizer] Chain {chain}: residue numbers start from 1.")
    # Note about only_pure_rna:
    # The flag and function `is_pure_rna` are kept for compatibility but intentionally not enforced.
    # This is to avoid unexpected early returns in diverse real-world structures.
    # In future, we may re-enable the filtering with a more robust detection.
    if only_pure_rna:
        # if not is_pure_rna(pdb_object, chain):
        #     print("The structure contains non-standard RNA bases or other molecules.")
        #     print("If you want to analyze them, please set pure_rna=False.")
        #     return
        print("[info] only_pure_rna flag is currently ignored (kept for compatibility).")
    
    # パーサーの選択に応じてベースペアを抽出
    if annotator.upper() == "DSSR":
        raw_df = dssr_wrapper(pdb_object, chain)
    elif annotator.upper() == "RNAVIEW":
        raw_df = rnaview_wrapper(pdb_object, chain)
    else:
        raise ValueError(f"Unsupported annotator: {annotator}. Use 'DSSR' or 'RNAView'.")
    
    processed_df = raw_df_processing(raw_df, annotator)
    processed_df, abnormal_pairs, dup_canonical_pairs = filter_abnormal_pairs(processed_df)
    
    # include_all フラグに基づいてフィルタリング
    if include_all:
        # すべての塩基対を使用（フィルタリング済み）
        filtered_df = processed_df
        canonical_count = len(processed_df[processed_df["is_canonical"]])
        noncanonical_count = len(processed_df) - canonical_count
        print(f"Using all base pairs: {len(filtered_df)} total ({canonical_count} canonical, {noncanonical_count} non-canonical)")
        # 重複する canonical pairs がある場合は警告
        if dup_canonical_pairs:
            print(f"Warning: {len(dup_canonical_pairs)} duplicate canonical pairs found and recorded.")
    else:
        # canonical base pairsのみを使用
        filtered_df = processed_df[processed_df["is_canonical"]]
        print(f"Using canonical base pairs only: {len(filtered_df)}/{len(processed_df)} pairs")
    
    # print(processed_df)
    # Orientation normalization: ensure i < j for PKextractor compatibility
    BPL = []
    for _, row in filtered_df.iterrows():
        i, j = row["position"]
        BPL.append((i, j) if i < j else (j, i))
    # print(f"extracted base pairs: {BPL}")
    PKlayers = PKextractor(BPL)

    if not skip_precoloring:
        print(f"Precoloring all atoms to white since skip_precoloring is {skip_precoloring}")
        cmd.color("white", f"{pdb_object} and chain {chain}")
    
    # Details dict for base pair analysis (similar to analysis script)
    records = filtered_df.to_dict(orient="records") if not filtered_df.empty else []
    details_dict = {tuple(rec["position"]): {
            "position": rec["position"],
            "residues": rec["residues"],
            "is_canonical": rec["is_canonical"],
            "saenger_id": rec["saenger_id"],
        } for rec in records
    }
    
    for depth, PKlayer in enumerate(PKlayers):
        color = get_color_for_depth(depth + 1, colors)
        
        # Layer statistics (when using non-canonical pairs)
        if include_all:
            canon_count = sum(1 for bp in PKlayer if details_dict[bp]["is_canonical"])
            noncanon_count = sum(1 for bp in PKlayer if not details_dict[bp]["is_canonical"])
            print(f"Layer {depth + 1}: {len(PKlayer)} pairs ({canon_count} canonical, {noncanon_count} non-canonical)")
        else:
            print(f"Layer {depth + 1}: {len(PKlayer)} canonical pairs")
        
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
            # Create selections with paper-friendly names and legacy names
            # paper-friendly: core, pk1, pk2, ... (core corresponds to layer 1)
            paper_name = "core" if depth == 0 else f"pk{depth}"
            paper_name = f"{str(pathlib.Path(pdb_object).stem)}_{chain}_{paper_name}"

            cmd.select(paper_name, f"{pdb_object} and chain {chain} and resi {selection_str}")
            print(f"Created selection: {paper_name} with residues {selection_str}")

            # legacy name kept for backward-compatibility
            # legacy_name = f"{pdb_object}_c{chain}_l{depth + 1}"
            # cmd.select(legacy_name, f"{pdb_object} and chain {chain} and resi {selection_str}")
            # print(f"Created selection: {legacy_name} with residues {selection_str}")
        print(f"Layer {depth + 1}: (i, j) = {PKlayer}")
    print("Coloring done.")
    print(f"pseudoknot order (number of layers): {len(PKlayers)}")
    
    clear_intermediate_files()
    return


cmd.extend("PseudoKnotVisualizer", PseudoKnotVisualizer)
cmd.extend("pkv", PseudoKnotVisualizer)

if __name__ == "__main__":
    pass
