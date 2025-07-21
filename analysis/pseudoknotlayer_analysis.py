"""
Pseudoknot Layer Analysis Script

このスクリプトはBGSUデータセットの全PDBファイルに対して：
1. RNAView or DSSRで塩基対情報を取得
2. PKextractorを使ってpseudoknot layerに分解
3. 各レイヤーでのcanonical/non-canonical base pairの割合を計算

$ cd PseudoknotVisualizer
$ python analysis/pseudoknotlayer_analysis.py --parser [RNAView or DSSR] [--canonical-only]

Date: 2025年7月21日
"""

import sys
import json
from pathlib import Path
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
    
    Args:
        pdb_file (Path): PDBファイルのパス
        parser (str): "RNAView" or "DSSR"
        
    Returns:
        dict: 解析結果
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
    all_bp_details_filtered, _, self_pairs_all = filter_self_pairs(all_bp_details, [(bp["position"][0], bp["position"][1]) for bp in all_bp_details])
    canonical_bp_details_filtered, canonical_bp_list_filtered, self_pairs = filter_self_pairs(canonical_bp_details, canonical_bp_list)
    
    has_self_pairs = len(self_pairs) > 0
    
    if has_self_pairs: 
        print(f"Warning: Found {len(self_pairs)} self-pairs (i=j) in {pdb_file.name}: {self_pairs}")
    
    if not all_bp_details_filtered:
        print(f"No valid base pairs found for {pdb_file.name} after removing self-pairs")
    
    # PKextractorでレイヤー分解（canonical BPのみ使用、自己ペアを除外）
    print("Running PKextractor...")

    pk_layers = PKextractor(canonical_bp_list_filtered.copy())

    print(f"Total pseudoknot layers extracted: {len(pk_layers)}")
    
    # 各レイヤーの解析
    layer_analysis = []

    bp_position_dict = {
        tuple(bp_detail["position"]): bp_detail 
        for bp_detail in canonical_bp_details_filtered
    }
    
    for layer_id, layer_bps in enumerate(pk_layers):
        # このレイヤーに含まれる塩基対の詳細情報を取得（自己ペアを除外）
        layer_bp_details = [
            bp_position_dict[bp_pos] for bp_pos in layer_bps if bp_pos in bp_position_dict
        ]

        # このレイヤーの塩基対は全てcanonical（PKextractorがcanonical BPのみで動作するため）
        layer_canonical_count = len(layer_bp_details)
        
        layer_analysis.append({
            "layer_id": layer_id,
            "total_bp_count": layer_canonical_count,
            "canonical_bp_count": layer_canonical_count,
            "non_canonical_bp_count": 0,
            # "base_pairs": layer_bp_details
        })
        
    result = {
        "pdb_id": pdb_file.stem,
        "chain_id": display_chain_id,  # ファイル名から抽出したチェーンID（表示用）
        "actual_chain_id": actual_chain_id,  # 実際に処理で使用したチェーンID
        "parser": parser,
        "total_bp_count": len(all_bp_details_filtered),
        "total_canonical_bp_count": len(canonical_bp_details_filtered),
        "self_pairs_count": len(self_pairs),
        "pseudoknot_layer_count": len(pk_layers),
        "layers": layer_analysis,
        "all_base_pairs": all_bp_details_filtered,  # all pairs
    }
    return result

def main():
    """メイン実行関数"""
    # コマンドライン引数を解析
    args = parse_args()
    
    # PDBファイルリストを取得
    pdb_files = get_pdb_files(DATASET_DIR) # Path object の list
    if not pdb_files:
        raise ValueError("No PDB files found in the dataset directory.")

    # 結果を格納するリスト
    results = []

    for pdb_file in tqdm(pdb_files, desc="Processing PDB files", unit="file"):
        result = analyze_single_pdb(pdb_file, args.parser, args.canonical_only)
        # 解析結果がNoneの場合はエラーを投げる
        if result is None:
            raise ValueError(f"Failed to analyze {pdb_file.name}")
        results.append(result)
    
    # 結果をJSONファイルに保存
    output_file = f"analysis/pseudoknot_analysis_{args.parser.lower()}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
