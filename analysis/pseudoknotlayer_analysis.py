"""
Pseudoknot Layer Analysis Script (Parallelized)

このスクリプトはBGSUデータセットの全PDBファイルに対して：
1. RNAView or DSSRで塩基対情報を取得
2. PKextractorを使ってpseudoknot layerに分解
3. 各レイヤーでのcanonical / non-canonical base pairの割合を計算

$ cd PseudoknotVisualizer
$ python analysis/pseudoknotlayer_analysis.py --parser [RNAView or DSSR] [--canonical-only]
"""

import sys
import json
from pathlib import Path
from functools import partial
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

script_dir = Path(__file__).parent.parent
sys.path.insert(0, str(script_dir))

# 分割したモジュールをインポート
from analysis.io_utils import (
    get_pdb_files,
    extract_chain_from_filename,
    extract_actual_chain_from_pdb,
    run_parser_analysis
)
from analysis.parsers import (
    parse_output_file,
    filter_self_pairs
)
from analysis.argparser import parse_args
from rna import PKextractor

# データセットディレクトリ
DATASET_DIR = "analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396"

def analyze_single_pdb(pdb_file, parser="RNAView", canonical_only=True):
    """
    単一のPDBファイルを解析
    """
    # チェーン情報を抽出（表示用）
    display_chain_id = extract_chain_from_filename(pdb_file.name)
    # PDBファイルのREMARK 350から実際のチェーンIDを取得
    actual_chain_id = extract_actual_chain_from_pdb(pdb_file)
    # パーサーを実行して出力ファイルを取得
    output_file = run_parser_analysis(pdb_file, actual_chain_id, parser)

    if output_file is None:
        print(f"Warning: {parser} output not found for {pdb_file.name}")
        raise ValueError(f"{parser} output not found for {pdb_file.name}")

    # 出力ファイルを解析して共通フォーマットで取得
    all_bp_details, canonical_bp_details, canonical_bp_list = parse_output_file(output_file, parser)

    # 自己ペア（i=j）を検出・除外
    all_bp_filtered, _, self_pairs_all = filter_self_pairs(
        all_bp_details,
        [(bp["position"][0], bp["position"][1]) for bp in all_bp_details]
    )
    can_bp_filtered, can_bp_list_filtered, self_pairs = filter_self_pairs(
        canonical_bp_details,
        canonical_bp_list
    )

    # レイヤー分解
    if canonical_only:
        pk_layers = PKextractor(can_bp_list_filtered.copy())
        bp_pos_dict = {tuple(bp["position"]): bp for bp in can_bp_filtered}
    else:
        print("Hello")
        pk_layers = PKextractor([(bp["position"][0], bp["position"][1]) for bp in all_bp_filtered.copy()])
        bp_pos_dict = {tuple(bp["position"]): bp for bp in all_bp_filtered}
    print("layer decomposed")

    layer_analysis = []
    for layer_id, layer_bps in enumerate(pk_layers):
        details = [bp_pos_dict[pos] for pos in layer_bps if pos in bp_pos_dict]
        if canonical_only:
            canon_count = len(layer_bps)
            noncanon_count = 0
        else:
            canon_count = sum(1 for bp in details if bp["is_canonical"])
            noncanon_count = sum(1 for bp in details if not bp["is_canonical"])
        layer_analysis.append({
            "layer_id": layer_id,
            "total_bp_count": len(details),
            "canonical_bp_count": canon_count,
            "non_canonical_bp_count": noncanon_count,
            "base_pairs": details
        })

    return {
        "pdb_id": pdb_file.stem,
        "chain_id": display_chain_id,
        "actual_chain_id": actual_chain_id,
        "parser": parser,
        "total_bp_count": len(all_bp_filtered),
        "total_canonical_bp_count": len(can_bp_filtered),
        "self_pairs_count": len(self_pairs),
        "pseudoknot_layer_count": len(pk_layers),
        "layers": layer_analysis,
        "all_base_pairs": all_bp_filtered # filtered: removed self-pairs
    }


def main():
    args = parse_args()
    pdb_files = get_pdb_files(DATASET_DIR)
    pdb_files = [Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396/1O9M_1_A-B.pdb")]

    if not pdb_files:
        raise ValueError("No PDB files found in the dataset directory.")

    # プロセス数指定がなければCPUコア数を使用
    n_procs = args.processes if hasattr(args, 'processes') and args.processes else cpu_count()
    n_procs = min(n_procs, 1)  

    # 並列で解析
    process_func = partial(
        analyze_single_pdb,
        parser=args.parser,
        canonical_only=args.canonical_only
    )
    process_func(pdb_files[0])  # テスト用に最初のファイルだけ実行
    return
    with Pool(processes=n_procs) as pool:
        results = list(tqdm(
            pool.imap(process_func, pdb_files),
            total=len(pdb_files),
            desc="Processing PDB files",
            unit="file"
        ))

    # 結果をJSONに保存
    output_file = f"analysis/pseudoknot_analysis_{args.parser.lower()}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()
