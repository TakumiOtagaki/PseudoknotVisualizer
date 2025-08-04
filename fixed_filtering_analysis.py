#!/usr/bin/env python3
"""
Fixed Detailed Filtering Analysis Script

RNAViewとDSSRの解析結果を正しく比較分析するスクリプト
"""

import json
from pathlib import Path
from collections import Counter

def analyze_parser_results(json_file, description, initial_chain_count):
    """各パーサーの結果を分析"""
    
    if not Path(json_file).exists():
        print(f"❌ {description}: ファイルが見つかりません ({json_file})")
        return None
    
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
    
    # Step 3: 解析成功の判定
    # output_existsフィールドがある場合（DSSR）とない場合（RNAView）で判定方法を変える
    if results and 'output_exists' in results[0]:
        # DSSRの場合
        successful_analyses = [r for r in results if r.get('output_exists', False)]
        success_criteria = "output_exists=True"
    else:
        # RNAViewの場合：total_bp_countが定義されているものを成功とする
        successful_analyses = [r for r in results if 'total_bp_count' in r]
        success_criteria = "total_bp_count フィールド存在"
    
    analysis_failures = processed_count - len(successful_analyses)
    print(f"Step 3 - 解析成功: {len(successful_analyses)} (失敗: {analysis_failures})")
    print(f"         判定基準: {success_criteria}")
    
    # Step 4: 塩基対が見つかったファイル
    chains_with_basepairs = [r for r in results if r.get('total_bp_count', 0) > 0]
    no_basepairs = len(successful_analyses) - len(chains_with_basepairs)
    print(f"Step 4 - 塩基対発見: {len(chains_with_basepairs)} (塩基対なし: {no_basepairs})")
    
    # 基本統計
    analyze_basic_stats(results, chains_with_basepairs, description)
    
    # 保持率計算
    retention_rate = len(chains_with_basepairs) / initial_chain_count * 100
    print(f"最終保持率: {retention_rate:.1f}%")
    print()
    
    return {
        'processed_count': processed_count,
        'successful_count': len(successful_analyses),
        'basepair_count': len(chains_with_basepairs),
        'results': results,
        'chains_with_basepairs': chains_with_basepairs
    }

def analyze_basic_stats(results, chains_with_basepairs, description):
    """基本統計の分析"""
    
    # 異常ペア統計
    if results:
        chains_with_abnormal = len([r for r in results if r.get('abnormal_pairs', [])])
        chains_with_duplicates = len([r for r in results if r.get('dup_canonical_pairs', [])])
        
        total_abnormal = sum(len(r.get('abnormal_pairs', [])) for r in results)
        total_duplicates = sum(len(r.get('dup_canonical_pairs', [])) for r in results)
        
        print(f"異常な塩基対が検出されたチェーン数: {chains_with_abnormal}")
        print(f"重複canonical塩基対が検出されたチェーン数: {chains_with_duplicates}")
        print(f"総異常ペア数: {total_abnormal}")
        print(f"総重複canonicalペア数: {total_duplicates}")
    
    # 塩基対数統計
    if chains_with_basepairs:
        bp_counts = [r['total_bp_count'] for r in chains_with_basepairs]
        canonical_bp_counts = [r.get('total_canonical_bp_count', 0) for r in chains_with_basepairs]
        
        print(f"平均塩基対数: {sum(bp_counts) / len(bp_counts):.1f}")
        print(f"平均canonical塩基対数: {sum(canonical_bp_counts) / len(canonical_bp_counts):.1f}")
        print(f"最大塩基対数: {max(bp_counts)}")
        print(f"最小塩基対数: {min(bp_counts)}")

def detailed_filtering_analysis():
    """フィルタリングステップの詳細分析"""
    
    # 初期データ
    dataset_dir = Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396")
    pdb_files = list(dataset_dir.glob("*.pdb"))
    initial_chain_count = len(pdb_files)
    
    print("=" * 80)
    print("フィルタリング詳細分析 (修正版)")
    print("=" * 80)
    print(f"初期PDBファイル数: {initial_chain_count}")
    print()
    
    # 解析結果ファイルを分析
    analysis_files = [
        ("analysis/pseudoknot_analysis_dssr_all.json", "DSSR (All)"),
        ("analysis/pseudoknot_analysis_rnaview_all.json", "RNAView (All)"),
        ("analysis/pseudoknot_analysis_dssr_canonical_only.json", "DSSR (Canonical Only)"),
        ("analysis/pseudoknot_analysis_rnaview_canonical_only.json", "RNAView (Canonical Only)")
    ]
    
    results_summary = {}
    
    for json_file, description in analysis_files:
        result = analyze_parser_results(json_file, description, initial_chain_count)
        if result:
            results_summary[description] = result
    
    # 比較分析
    print("🔍 パーサー比較分析")
    print("=" * 80)
    
    if "DSSR (All)" in results_summary and "RNAView (All)" in results_summary:
        dssr_result = results_summary["DSSR (All)"]
        rnaview_result = results_summary["RNAView (All)"]
        
        print("DSSRとRNAViewの比較:")
        print(f"DSSR   - 処理成功: {dssr_result['successful_count']}, 塩基対発見: {dssr_result['basepair_count']}")
        print(f"RNAView - 処理成功: {rnaview_result['successful_count']}, 塩基対発見: {rnaview_result['basepair_count']}")
        
        # なぜRNAViewの解析成功数が0だったかの説明
        print("\n💡 RNAViewで解析成功数が0だった理由:")
        print("- DSSRの結果にはoutput_existsフィールドがあり、これで解析成功を判定")
        print("- RNAViewの結果にはoutput_existsフィールドがない")
        print("- 修正後は、total_bp_countフィールドの存在で解析成功を判定")

if __name__ == "__main__":
    detailed_filtering_analysis()
