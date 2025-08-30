#!/usr/bin/env python3
"""
Unified Filtering Analysis Script

BGSU__M__All__A__4_0__pdb_3_396データセットの解析結果を統計的にまとめる統合スクリプト
"""

import json
from pathlib import Path
import pandas as pd

def analyze_parser_results(json_file, description, initial_count):
    """パーサー結果の分析"""
    with open(json_file, 'r') as f:
        results = json.load(f)
    
    processed_count = len(results)
    
    # 解析成功の判定（DSSRはoutput_exists、RNAViewはtotal_bp_countフィールドの存在）
    if results and 'output_exists' in results[0]:
        successful = [r for r in results if r.get('output_exists', False)]
    else:
        successful = [r for r in results if 'total_bp_count' in r]
    
    chains_with_bp = [r for r in results if r.get('total_bp_count', 0) > 0]
    
    # 統計計算
    if chains_with_bp:
        bp_counts = [r['total_bp_count'] for r in chains_with_bp]
        canonical_counts = [r.get('total_canonical_bp_count', 0) for r in chains_with_bp]
        avg_bp = sum(bp_counts) / len(bp_counts)
        avg_canonical = sum(canonical_counts) / len(canonical_counts)
        max_bp = max(bp_counts)
        min_bp = min(bp_counts)
    else:
        avg_bp = avg_canonical = max_bp = min_bp = 0
    
    # 異常ペア統計
    abnormal_chains = len([r for r in results if r.get('abnormal_pairs', [])])
    chains_with_dup_canonical = len([r for r in results if r.get('dup_canonical_pairs', [])])
    total_abnormal = sum(len(r.get('abnormal_pairs', [])) for r in results)
    total_dup_canonical = sum(len(r.get('dup_canonical_pairs', [])) for r in results)
    
    retention_rate = len(chains_with_bp) / initial_count * 100
    
    return {
        'description': description,
        'processed': processed_count,
        'successful': len(successful),
        'with_basepairs': len(chains_with_bp),
        'abnormal_chains': abnormal_chains,
        'chains_with_dup_canonical': chains_with_dup_canonical,
        'total_abnormal': total_abnormal,
        'total_dup_canonical': total_dup_canonical,
        'avg_bp': avg_bp,
        'avg_canonical': avg_canonical,
        'max_bp': max_bp,
        'min_bp': min_bp,
        'retention_rate': retention_rate
    }

def main():
    """メイン分析関数"""
    # 初期データ
    dataset_dir = Path("analysis/datasets/BGSU__M__All__A__4_0__pdb_3_396")
    initial_count = len(list(dataset_dir.glob("*.pdb")))
    
    print("=" * 80)
    print("BGSU__M__All__A__4_0__pdb_3_396 統合フィルタリング分析")
    print("=" * 80)
    print(f"初期PDBファイル数: {initial_count}")
    print()
    
    # 解析対象ファイル
    analysis_files = [
        ("analysis/pseudoknot_analysis_dssr_all.json", "DSSR (All)"),
        ("analysis/pseudoknot_analysis_dssr_canonical_only.json", "DSSR (Canonical Only)"),
        ("analysis/pseudoknot_analysis_rnaview_all.json", "RNAView (All)"),
        ("analysis/pseudoknot_analysis_rnaview_canonical_only.json", "RNAView (Canonical Only)")
    ]
    
    results = []
    for json_file, description in analysis_files:
        if Path(json_file).exists():
            result = analyze_parser_results(json_file, description, initial_count)
            results.append(result)
    
    # 結果表示
    df_data = []
    for r in results:
        df_data.append({
            'parser': r['description'],
            'processed': r['processed'],
            'successful': r['successful'],
            'with_basepairs': r['with_basepairs'],
            'abnormal_chains': r['abnormal_chains'],
            'dup_canonical_chains': r['chains_with_dup_canonical'],
            'retention_rate(%)': f"{r['retention_rate']:.1f}"
        })
    
    df = pd.DataFrame(df_data)
    print("📊 フィルタリング結果サマリー")
    print(df.to_string(index=False))
    print()
    
    # 詳細統計
    for r in results:
        print(f"📊 {r['description']}")
        print(f"   フロー: {initial_count} → {r['processed']} → {r['successful']} → {r['with_basepairs']}")
        print(f"   塩基対: 平均{r['avg_bp']:.1f} (canonical{r['avg_canonical']:.1f}) 範囲[{r['min_bp']}-{r['max_bp']}]")
        print(f"   異常ペア: {r['total_abnormal']} 重複canonical: {r['total_dup_canonical']}")
        print()
    
    # フィルタリング説明
    print("🔍 実装フィルタリング")
    print("-" * 50)
    print("1. 事前除外: 問題ファイル (PDB_00003OK4_1_2.pdb)")
    print("2. 解析失敗: DSSR/RNAView実行エラー")
    print("3. 塩基対なし: 構造内に塩基対が検出されない")
    print("4. 異常ペア除外:")
    print("   - 自己ペア (i=i)")
    print("   - 方向重複 (i,j)と(j,i)")
    print("   - 位置重複 (i,j)と(i,j') → canonicalを優先")
    print("5. Canonical-only: canonical塩基対のみ使用")
    print()
    print("💡 dup_canonical_chains の説明:")
    print("   同一チェーン内で同じ塩基位置に複数のcanonical塩基対が")
    print("   検出されたチェーンの数（例: 位置5-15に2つのcanonical塩基対）")

if __name__ == "__main__":
    main()
