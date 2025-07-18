"""
Pseudoknot Layer Analysis Script

このスクリプトはBGSUデータセットの全PDBファイルに対して：
1. RNAView or DSSRで塩基対情報を取得
2. PKextractorを使ってpseudoknot layerに分解
3. 各レイヤーでのcanonical/non-canonical base pairの割合を計算

Author: PseudoknotVisualizer Analysis
Date: 2025年7月18日
"""

import os
import sys
import glob
import json
import traceback
from pathlib import Path
from tqdm import tqdm

script_dir = Path(__file__).parent.parent
print(f"Script directory: {script_dir}")
sys.path.insert(0, str(script_dir))

from addressRNAviewOutput import load_rnaview_data, canonical_extraction_from_rnaview_df
from addressDSSROutput import load_dssr_data, canonical_extraction_from_dssr_df
from CLI_PseudoknotVisualizer import CLI_rnaview, CLI_dssr
from rna import PKextractor
import pandas as pd

# データセットディレクトリ
DATASET_DIR = "analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396"

def get_pdb_files():
    """データセットディレクトリから全PDBファイルのリストを取得"""
    pdb_files = list(Path(DATASET_DIR).glob("*.pdb"))
    print(f"Found {len(pdb_files)} PDB files in dataset")
    return sorted(pdb_files)

def extract_chain_from_filename(filename):
    """
    ファイル名からチェーン情報を抽出
    例: 1EHZ_1_A.pdb -> A
    例: 1J7T_1_A-B.pdb -> A (最初のチェーンを使用)
    """
    stem = Path(filename).stem
    parts = stem.split('_')
    if len(parts) >= 3:
        chain_part = parts[2]
        # ハイフンで区切られている場合は最初のチェーンを使用
        chain = chain_part.split('-')[0]
        return chain
    return 'A'  # デフォルト

def analyze_single_pdb(pdb_file, parser="RNAView"):
    """
    単一のPDBファイルを解析
    
    Args:
        pdb_file (Path): PDBファイルのパス
        parser (str): "RNAView" or "DSSR"
        
    Returns:
        dict: 解析結果
    """
    try:
        # チェーン情報を抽出
        chain_id = extract_chain_from_filename(pdb_file.name)
        
        print(f"Processing {pdb_file.name} with chain {chain_id} using {parser}...")
        
        # パーサーに応じてベースペア情報を取得
        if parser.upper() == "RNAVIEW":
            # RNAViewで全ベースペア情報を取得
            try:
                CLI_rnaview(str(pdb_file), chain_id)
                # 出力ファイルパスを構築
                intermediate_dir = parent_dir / "intermediate"
                output_file = intermediate_dir / f"{pdb_file.name}.out"
                
                if not output_file.exists():
                    print(f"Warning: RNAView output not found for {pdb_file.name}")
                    return None
                    
                # 全データをロード
                all_bp_df = load_rnaview_data(str(output_file))
                canonical_bp_df = canonical_extraction_from_rnaview_df(all_bp_df)
                
            except Exception as e:
                print(f"Error processing {pdb_file.name} with RNAView: {e}")
                return None
                
        elif parser.upper() == "DSSR":
            # DSSRで全ベースペア情報を取得
            try:
                CLI_dssr(str(pdb_file), chain_id)
                # 出力ファイルパスを構築
                intermediate_dir = parent_dir / "intermediate"
                output_file = intermediate_dir / f"{pdb_file.name}.dssr.json"
                
                if not output_file.exists():
                    print(f"Warning: DSSR output not found for {pdb_file.name}")
                    return None
                
                # 全データをロード
                all_bp_df = load_dssr_data(str(output_file))
                canonical_bp_df = canonical_extraction_from_dssr_df(all_bp_df)
                
            except Exception as e:
                print(f"Error processing {pdb_file.name} with DSSR: {e}")
                return None
        else:
            raise ValueError(f"Unsupported parser: {parser}")
        
        # ベースペアリストの準備
        all_bp_list = [(row["left_idx"], row["right_idx"]) for _, row in all_bp_df.iterrows()]
        canonical_bp_list = [(row["left_idx"], row["right_idx"]) for _, row in canonical_bp_df.iterrows()]
        
        if not all_bp_list:
            print(f"No base pairs found for {pdb_file.name}")
            return None
        
        # PKextractorでレイヤー分解（canonical BPのみ使用）
        pk_layers = PKextractor(canonical_bp_list.copy()) if canonical_bp_list else []
        
        # 各レイヤーの解析
        layer_analysis = []
        total_canonical = len(canonical_bp_list)
        total_all = len(all_bp_list)
        total_non_canonical = total_all - total_canonical
        
        for layer_id, layer_bps in enumerate(pk_layers):
            layer_canonical_count = len(layer_bps)
            # このレイヤーの塩基対は全てcanonical（PKextractorがcanonical BPのみで動作するため）
            layer_non_canonical_count = 0
            
            layer_analysis.append({
                "layer_id": layer_id,
                "total_bp_count": layer_canonical_count,
                "canonical_bp_count": layer_canonical_count,
                "non_canonical_bp_count": layer_non_canonical_count,
                "canonical_ratio": 1.0,  # 100%
                "non_canonical_ratio": 0.0,
                "base_pairs": layer_bps
            })
        
        result = {
            "pdb_id": pdb_file.stem,
            "chain_id": chain_id,
            "parser": parser,
            "total_bp_count": total_all,
            "total_canonical_bp_count": total_canonical,
            "total_non_canonical_bp_count": total_non_canonical,
            "total_canonical_ratio": total_canonical / total_all if total_all > 0 else 0,
            "total_non_canonical_ratio": total_non_canonical / total_all if total_all > 0 else 0,
            "pseudoknot_layer_count": len(pk_layers),
            "layers": layer_analysis
        }
        
        return result
        
    except Exception as e:
        print(f"Error analyzing {pdb_file.name}: {e}")
        print(traceback.format_exc())
        return None

def main():
    """メイン実行関数"""
    print("=== Pseudoknot Layer Analysis ===")
    print(f"Dataset directory: {DATASET_DIR}")
    
    # PDBファイルリストを取得
    pdb_files = get_pdb_files()
    print(f"Found {len(pdb_files)} PDB files to analyze.")
    
    if not pdb_files:
        print("No PDB files found!")
        return
    
    # 解析設定
    parser = "RNAView"  # "RNAView" or "DSSR"
    
    print(f"Using parser: {parser}")
    print(f"Processing {len(pdb_files)} files...")
    
    # 結果を格納するリスト
    results = []
    successful_count = 0
    failed_count = 0
    
    # 各PDBファイルを処理（最初の10個をテスト用に制限）
    # for i, pdb_file in enumerate(pdb_files[:10]):  # テスト用に最初の10個のみ
    for i, pdb_file in enumerate(tqdm(pdb_files, desc="Processing PDB files", unit="file")):
        print(f"\n--- Processing {i+1}/{min(10, len(pdb_files))}: {pdb_file.name} ---")
        
        result = analyze_single_pdb(pdb_file, parser)
        
        if result is not None:
            results.append(result)
            successful_count += 1
            print(f"Successfully processed {pdb_file.name}")
            print(f"Total BP: {result['total_bp_count']}, Canonical: {result['total_canonical_bp_count']}, Layers: {result['pseudoknot_layer_count']}")
        else:
            failed_count += 1
            print(f"Failed to process {pdb_file.name}")
    
    # 結果をJSONファイルに保存
    output_file = f"analysis/pseudoknot_analysis_{parser.lower()}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print(f"\n=== Analysis Complete ===")
    print(f"Successfully processed: {successful_count}")
    print(f"Failed: {failed_count}")
    print(f"Results saved to: {output_file}")
    
    # 簡単な統計を表示
    if results:
        total_structures = len(results)
        avg_total_bp = sum(r['total_bp_count'] for r in results) / total_structures
        avg_canonical_ratio = sum(r['total_canonical_ratio'] for r in results) / total_structures
        avg_layer_count = sum(r['pseudoknot_layer_count'] for r in results) / total_structures
        
        print(f"\n=== Summary Statistics ===")
        print(f"Average total base pairs per structure: {avg_total_bp:.1f}")
        print(f"Average canonical ratio: {avg_canonical_ratio:.3f}")
        print(f"Average pseudoknot layer count: {avg_layer_count:.1f}")

if __name__ == "__main__":
    main()
