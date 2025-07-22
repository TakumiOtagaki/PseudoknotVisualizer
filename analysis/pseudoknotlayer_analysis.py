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
import pandas as pd
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
    filter_abnormal_pairs,
    raw_df_processing
)
from analysis.argparser import parse_args
from rna import PKextractor

# データセットディレクトリ
DATASET_DIR = "analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396"
remove_files = [
    "/large/otgk/PseudoknotVisualizer/intermediate/PDB_00003OK4_1_2.pdb"
]

def analyze_single_pdb(pdb_file, parser="RNAView", canonical_only=True):
    """
    単一のPDBファイルを解析
    """
    # チェーン情報を抽出（表示用）
    display_chain_id = extract_chain_from_filename(pdb_file.name)
    # PDBファイルのREMARK 350から実際のチェーンIDを取得
    actual_chain_id = extract_actual_chain_from_pdb(pdb_file)
    # パーサーを実行して出力ファイルを取得
    output_file, raw_df = run_parser_analysis(pdb_file, actual_chain_id, parser)
    print(f"Output file for {pdb_file.name}: {output_file}")

    if output_file is None:
        print(f"Warning: {parser} output not found for {pdb_file.name}")
        raise ValueError(f"{parser} output not found for {pdb_file.name}")

    # 出力ファイルを解析して共通フォーマットで取得
    # raw_df = parse_output_file(output_file, parser)
    print("raw_df:\n", raw_df.head())
    processed_df = raw_df_processing(raw_df, parser)
    print(f"Processed DataFrame for {pdb_file.name}:\n", processed_df)

    # 自己ペア（i=j）を検出・除外
    processed_df, abnormal_pairs, dup_canonical_pairs = filter_abnormal_pairs(processed_df)

    print(f"Analyzing {pdb_file.name} (Chain: {display_chain_id}, Actual Chain: {actual_chain_id})")
    print("processed_df:\n", processed_df.head())
    print(f"all base pairs found: {len(processed_df)}")
    # 以下のコードは、dict が空の時にエラーになる。defaultdict を使えば問題ないが、未実装です。

    canonical_processed_df = processed_df[processed_df["is_canonical"]] if not processed_df.empty else pd.DataFrame(columns=processed_df.columns)
    if canonical_only:
        # もし共通している (i, j) と (i, j') のような塩s基対があれば、error という扱いにして飛ばす
        basepair_list = [ (bp[0], bp[1]) for bp in canonical_processed_df["position"]]
    else:
        basepair_list = [ (bp[0], bp[1]) for bp in processed_df["position"]]
    # key: pair (i, j) , value: base pair details
    details_dict = {(bp[0][0], bp[0][1]): {
            "position": bp[0], "residues": bp, "is_canonical": bp[2], "saenger_id": bp[3]
        } for bp in processed_df.values.tolist()
    }
    # print(f"Total base pairs found: {details_dict}")
    print("base pair list:")
    for bp in basepair_list:
        print(f"  {bp} ")
    print("duplicated canonical pairs:")
    for bp in dup_canonical_pairs:
        print(f"  {bp} ")
    pk_layers = PKextractor(basepair_list.copy())
    print("layer decomposed")    

    layer_analysis = []
    for layer_id, layer_bps in enumerate(pk_layers):
        if canonical_only:
            canon_count = len(layer_bps)
            noncanon_count = 0
        else:
            canon_count = sum(1 for bp in layer_bps if details_dict[bp]["is_canonical"])
            noncanon_count = sum(1 for bp in layer_bps if not details_dict[bp]["is_canonical"])
        layer_analysis.append({
            "layer_id": layer_id,
            "total_bp_count": len(layer_bps),
            "basepair_details": [details_dict[bp] for bp in layer_bps],
            "canonical_bp_count": canon_count,
            "non_canonical_bp_count": noncanon_count,
        })

    return {
        "pdb_id": pdb_file.stem,
        "chain_id": display_chain_id,
        "actual_chain_id": actual_chain_id,
        "parser": parser,
        "total_bp_count": len(processed_df),
        "total_canonical_bp_count": len(canonical_processed_df),
        "pseudoknot_layer_count": len(pk_layers),
        "layers": layer_analysis,
        "all_base_pairs": processed_df["position"].tolist(),  # filtered: removed self-pairs
        "abnormal_pairs": abnormal_pairs,
        # "details": details_dict.values,
        "dup_canonical_pairs": list(dup_canonical_pairs) if dup_canonical_pairs else [],
    }


def main():
    args = parse_args()
    pdb_files = get_pdb_files(DATASET_DIR)
    pdb_files = [pdb_file for pdb_file in pdb_files if pdb_file not in remove_files]
    # pdb_files = [Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396/1O9M_1_A-B.pdb")]
    # pdb_files = [Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396/PDB_00003OK4_1_2.pdb")]
    # pdb_files = pdb_files[:10]
    if not pdb_files:
        raise ValueError("No PDB files found in the dataset directory.")

    # プロセス数指定がなければCPUコア数を使用
    n_procs = args.processes if hasattr(args, 'processes') and args.processes else cpu_count()
    n_procs = min(n_procs, args.ncpus)  

    # 並列で解析
    process_func = partial(
        analyze_single_pdb,
        parser=args.parser,
        canonical_only=args.canonical_only
    )
    process_func(pdb_files[0])  # テスト用に最初のファイルだけ実行
    # return

    # with Pool(processes=n_procs) as pool:
    #     results = list(tqdm(
    #         pool.imap(process_func, pdb_files),
    #         total=len(pdb_files),
    #         desc="Processing PDB files",
    #         unit="file"
    #     ))
    results = []
    for pdb_file in tqdm(pdb_files, desc="Processing PDB files", unit="file"):
        result = process_func(pdb_file)
        results.append(result)
        print(f"Processed {pdb_file.name}: {result['total_bp_count']} base pairs found.")
        print("-" * 40)

    # 結果をJSONに保存
    output_file = f"analysis/pseudoknot_analysis_{args.parser.lower()}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()
