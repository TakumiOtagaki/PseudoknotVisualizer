#!/usr/bin/env python3
"""
Detailed Filtering Analysis Script

各フィルタリングステップでのチェーン数の変化を詳細に追跡
"""

import json
from pathlib import Path
from collections import Counter

def detailed_filtering_analysis():
    """フィルタリングステップの詳細分析"""
    
    # 初期データ
    dataset_dir = Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396")
    pdb_files = list(dataset_dir.glob("*.pdb"))
    initial_chain_count = len(pdb_files)
    
    print("=" * 80)
    print("フィルタリング詳細分析")
    print("=" * 80)
    
    # 解析結果ファイルをすべて分析
    analysis_files = [
        ("analysis/pseudoknot_analysis_dssr_all.json", "DSSR (All)"),
        ("analysis/pseudoknot_analysis_rnaview_all.json", "RNAView (All)")
    ]
    
    for json_file, description in analysis_files:
        if not Path(json_file).exists():
            print(f"❌ {json_file} が見つかりません")
            continue
        
        with open(json_file, 'r') as f:
            results = json.load(f)
        
        print(f"📊 フィルタリングフロー分析 ({description})")
        print("-" * 80)
    
        # Step 1: 初期ファイル数
        print(f"Step 1 - 初期PDBファイル数: {initial_chain_count}")
        
        # Step 2: 処理対象ファイル数
        processed_count = len(results)
        removed_files = initial_chain_count - processed_count
        print(f"Step 2 - 処理対象ファイル数: {processed_count} (除外: {removed_files})")
        
        # Step 3: 解析成功 (output_existsフィールドがない場合は、total_bp_count >= 0で判定)
        if 'output_exists' in results[0] if results else False:
            successful_analyses = [r for r in results if r.get('output_exists', False)]
        else:
            # RNAViewの場合：total_bp_countが定義されているものを成功とする
            successful_analyses = [r for r in results if 'total_bp_count' in r]
        
        analysis_failures = processed_count - len(successful_analyses)
        print(f"Step 3 - 解析成功: {len(successful_analyses)} (失敗: {analysis_failures})")
        
        # Step 4: 塩基対が見つかったファイル
        chains_with_basepairs = [r for r in results if r.get('total_bp_count', 0) > 0]
        no_basepairs = len(successful_analyses) - len(chains_with_basepairs)
        print(f"Step 4 - 塩基対発見: {len(chains_with_basepairs)} (塩基対なし: {no_basepairs})")
        
        # フィルタリング詳細統計
        analyze_filtering_details(results, description)
        print()

def analyze_filtering_details(results, description):
    """フィルタリング詳細統計の分析"""
    print("\n🔍 フィルタリング詳細統計")
    print("-" * 80)    # 異常ペアの種類別統計
    self_pairs_count = 0
    duplicate_pairs_count = 0
    other_abnormal_count = 0
    
    for result in results:
        abnormal_pairs = result.get('abnormal_pairs', [])
        dup_canonical_pairs = result.get('dup_canonical_pairs', [])
        
        # 自己ペアをカウント
        for pair in abnormal_pairs:
            if isinstance(pair, list) and len(pair) == 2 and pair[0] == pair[1]:
                self_pairs_count += 1
            else:
                other_abnormal_count += 1
        
        duplicate_pairs_count += len(dup_canonical_pairs)
    
    print(f"自己ペア (i,i) 検出数: {self_pairs_count}")
    print(f"重複canonical塩基対検出数: {duplicate_pairs_count}")
    print(f"その他異常ペア検出数: {other_abnormal_count}")

def analyze_distribution_stats(results, description):
    """分布統計の分析"""
    chains_with_basepairs = [r for r in results if r.get('total_bp_count', 0) > 0]
    
    if not chains_with_basepairs:
        print("塩基対が見つかったチェーンがありません")
        return
    
    # 塩基対数分布
    print("\n📈 塩基対数分布")
    print("-" * 80)
    
    bp_counts = [r['total_bp_count'] for r in chains_with_basepairs]
    canonical_bp_counts = [r.get('total_canonical_bp_count', 0) for r in chains_with_basepairs]
    
    # ヒストグラム的表示
    bp_ranges = [(0, 10), (11, 20), (21, 50), (51, 100), (101, 200), (201, 500), (500, float('inf'))]
    
    for min_bp, max_bp in bp_ranges:
        if max_bp == float('inf'):
            count = len([bp for bp in bp_counts if bp > min_bp])
            range_str = f"{min_bp}+"
        else:
            count = len([bp for bp in bp_counts if min_bp <= bp <= max_bp])
            range_str = f"{min_bp}-{max_bp}"
        
        percentage = count / len(bp_counts) * 100 if bp_counts else 0
        print(f"塩基対数 {range_str:>8}: {count:>4} チェーン ({percentage:>5.1f}%)")
    
    # Canonical ratio分析
    print("\n🎯 Canonical塩基対比率分析")
    print("-" * 80)
    
    canonical_ratios = []
    for result in chains_with_basepairs:
        total_bp = result['total_bp_count']
        canonical_bp = result.get('total_canonical_bp_count', 0)
        if total_bp > 0:
            ratio = canonical_bp / total_bp
            canonical_ratios.append(ratio)
    
    if canonical_ratios:
        avg_canonical_ratio = sum(canonical_ratios) / len(canonical_ratios)
        high_canonical = len([r for r in canonical_ratios if r >= 0.8])
        medium_canonical = len([r for r in canonical_ratios if 0.5 <= r < 0.8])
        low_canonical = len([r for r in canonical_ratios if r < 0.5])
        
        print(f"平均canonical比率: {avg_canonical_ratio:.3f}")
        print(f"高canonical比率 (≥80%): {high_canonical} チェーン ({high_canonical/len(canonical_ratios)*100:.1f}%)")
        print(f"中canonical比率 (50-79%): {medium_canonical} チェーン ({medium_canonical/len(canonical_ratios)*100:.1f}%)")
        print(f"低canonical比率 (<50%): {low_canonical} チェーン ({low_canonical/len(canonical_ratios)*100:.1f}%)")
    
    # 最終サマリー
    print("\n📋 最終フィルタリングサマリー")
    print("-" * 80)
    retention_rates = [
        ("初期ファイル数", initial_chain_count, 100.0),
        ("処理対象", processed_count, processed_count/initial_chain_count*100),
        ("解析成功", len(successful_analyses), len(successful_analyses)/initial_chain_count*100),
        ("塩基対発見", len(chains_with_basepairs), len(chains_with_basepairs)/initial_chain_count*100),
    ]
    
    for stage, count, percentage in retention_rates:
        print(f"{stage:>12}: {count:>4} チェーン ({percentage:>5.1f}%)")

if __name__ == "__main__":
    detailed_filtering_analysis()
