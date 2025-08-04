#!/usr/bin/env python3
"""
Analysis Summary Script

BGSU__M__All__A__4_0__pdb_3_396データセットの解析結果を統計的にまとめるスクリプト
フィルタリング前後のチェーン数の変化を追跡します。
"""

import json
from pathlib import Path

def load_analysis_results(json_file):
    """解析結果JSONファイルを読み込み"""
    with open(json_file, 'r') as f:
        return json.load(f)

def analyze_filtering_steps():
    """フィルタリングステップの詳細分析"""
    
    # データセットディレクトリのPDBファイル数（フィルタリング前の総チェーン数）
    dataset_dir = Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396")
    pdb_files = list(dataset_dir.glob("*.pdb"))
    initial_chain_count = len(pdb_files)
    
    print("=" * 60)
    print("BGSU__M__All__A__4_0__pdb_3_396 データセット分析")
    print("=" * 60)
    print(f"フィルタリング前の総チェーン数: {initial_chain_count}")
    print()
    
    # 解析結果ファイルを読み込み
    analysis_files = [
        ("analysis/pseudoknot_analysis_dssr_all.json", "DSSR (All)"),
        ("analysis/pseudoknot_analysis_dssr_canonical_only.json", "DSSR (Canonical Only)"),
        ("analysis/pseudoknot_analysis_rnaview_all.json", "RNAView (All)"),
        ("analysis/pseudoknot_analysis_rnaview_canonical_only.json", "RNAView (Canonical Only)")
    ]
    
    for json_file, description in analysis_files:
        if not Path(json_file).exists():
            print(f"⚠️  {description}: ファイルが見つかりません ({json_file})")
            continue
            
        results = load_analysis_results(json_file)
        
        print(f"📊 {description}")
        print("-" * 50)
        
        # 基本統計
        total_processed = len(results)
        successful_chains = [r for r in results if r.get('output_exists', False)]
        chains_with_basepairs = [r for r in results if r.get('total_bp_count', 0) > 0]
        
        print(f"処理されたチェーン数: {total_processed}")
        print(f"解析成功チェーン数: {len(successful_chains)}")
        print(f"塩基対が見つかったチェーン数: {len(chains_with_basepairs)}")
        
        # フィルタリング詳細
        chains_with_abnormal = [r for r in results if r.get('abnormal_pairs', [])]
        chains_with_duplicates = [r for r in results if r.get('dup_canonical_pairs', [])]
        
        print(f"異常な塩基対が検出されたチェーン数: {len(chains_with_abnormal)}")
        print(f"重複canonical塩基対が検出されたチェーン数: {len(chains_with_duplicates)}")
        
        # 塩基対統計
        if chains_with_basepairs:
            total_bp_counts = [r['total_bp_count'] for r in chains_with_basepairs]
            canonical_bp_counts = [r.get('total_canonical_bp_count', 0) for r in chains_with_basepairs]
            
            print(f"平均塩基対数: {sum(total_bp_counts) / len(total_bp_counts):.1f}")
            print(f"平均canonical塩基対数: {sum(canonical_bp_counts) / len(canonical_bp_counts):.1f}")
            print(f"最大塩基対数: {max(total_bp_counts)}")
            print(f"最小塩基対数: {min(total_bp_counts)}")
        
        print()
    
    # フィルタリングの詳細分析
    print("🔍 フィルタリング詳細分析")
    print("=" * 60)
    
    # DSSR Allの結果を使って詳細分析
    dssr_all_file = "analysis/pseudoknot_analysis_dssr_all.json"
    if Path(dssr_all_file).exists():
        results = load_analysis_results(dssr_all_file)
        
        print("DSSRで実装されているフィルタリング:")
        print("1. 自己ペア (i=i) の除外")
        print("2. 重複した位置の塩基対の処理:")
        print("   - (i,j) と (i,j') のような重複があった場合:")
        print("     • 一方がcanonical、他方がnon-canonicalなら、canonicalを保持")
        print("     • 両方ともnon-canonicalなら、両方とも除外")
        print("     • 両方ともcanonicalなら、重複canonicalペアとして記録")
        print("3. 方向の統一 (i,j) と (j,i) の重複処理")
        print()
        
        # 異常ペアの詳細統計
        total_abnormal_pairs = 0
        total_duplicate_canonical = 0
        
        for result in results:
            abnormal_count = len(result.get('abnormal_pairs', []))
            duplicate_count = len(result.get('dup_canonical_pairs', []))
            total_abnormal_pairs += abnormal_count
            total_duplicate_canonical += duplicate_count
        
        print(f"総異常ペア数: {total_abnormal_pairs}")
        print(f"総重複canonicalペア数: {total_duplicate_canonical}")
        
        # チェーン数の変化を追跡
        print("\n📈 フィルタリング後のチェーン数変化:")
        print(f"初期チェーン数: {initial_chain_count}")
        print(f"処理されたチェーン数: {len(results)}")
        print(f"解析成功チェーン数: {len([r for r in results if r.get('output_exists', False)])}")
        print(f"有効な塩基対があるチェーン数: {len([r for r in results if r.get('total_bp_count', 0) > 0])}")
        
        retention_rate = len([r for r in results if r.get('total_bp_count', 0) > 0]) / initial_chain_count * 100
        print(f"最終保持率: {retention_rate:.1f}%")

if __name__ == "__main__":
    analyze_filtering_steps()
