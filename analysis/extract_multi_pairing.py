#!/usr/bin/env python3
"""
Multi Base Pairing Entry Extractor

DSSR解析結果から"(i, j) と (i, j') といった multi base pairing を持ったエントリ"を抽出し、
新しいJSONファイルとして保存するスクリプト

Step by step で慎重に実行します。
"""

import json
from pathlib import Path

def extract_multi_pairing_entries(input_json_file, output_json_file):
    """
    multi base pairingを持つエントリを抽出
    
    Args:
        input_json_file (str): 入力JSONファイルパス
        output_json_file (str): 出力JSONファイルパス
    
    Returns:
        int: 抽出されたエントリ数
    """
    print(f"📖 入力ファイル読み込み: {input_json_file}")
    
    with open(input_json_file, 'r') as f:
        all_results = json.load(f)
    
    print(f"📊 総エントリ数: {len(all_results)}")
    
    # multi base pairingを持つエントリを抽出
    multi_pairing_entries = []
    
    for entry in all_results:
        dup_canonical_pairs = entry.get('dup_canonical_pairs', [])
        
        # dup_canonical_pairsが空でないエントリを抽出
        if dup_canonical_pairs:
            multi_pairing_entries.append(entry)
            print(f"✅ Multi pairing発見: {entry['pdb_id']} - {len(dup_canonical_pairs)} pairs")
    
    print(f"\n🎯 Multi base pairing エントリ数: {len(multi_pairing_entries)}")
    
    # 結果を保存
    with open(output_json_file, 'w') as f:
        json.dump(multi_pairing_entries, f, indent=2, ensure_ascii=False)
    
    print(f"💾 保存完了: {output_json_file}")
    
    return len(multi_pairing_entries)

def analyze_multi_pairing_stats(entries):
    """
    multi base pairingエントリの統計を分析
    """
    print("\n📈 Multi Base Pairing 統計分析")
    print("=" * 50)
    
    total_dup_pairs = 0
    layer_counts = []
    bp_counts = []
    canonical_ratios = []
    
    for entry in entries:
        dup_pairs = len(entry.get('dup_canonical_pairs', []))
        total_dup_pairs += dup_pairs
        
        layer_count = entry.get('pseudoknot_layer_count', 0)
        layer_counts.append(layer_count)
        
        total_bp = entry.get('total_bp_count', 0)
        canonical_bp = entry.get('total_canonical_bp_count', 0)
        bp_counts.append(total_bp)
        
        if total_bp > 0:
            canonical_ratios.append(canonical_bp / total_bp)
    
    print(f"総重複canonical塩基対数: {total_dup_pairs}")
    print(f"平均重複数/エントリ: {total_dup_pairs / len(entries):.1f}")
    print(f"平均pseudoknot layer数: {sum(layer_counts) / len(layer_counts):.1f}")
    print(f"平均塩基対数: {sum(bp_counts) / len(bp_counts):.1f}")
    print(f"平均canonical比率: {sum(canonical_ratios) / len(canonical_ratios):.3f}")
    
    # レイヤー数分布
    from collections import Counter
    layer_dist = Counter(layer_counts)
    print(f"\nPseudoknot Layer分布:")
    for layers, count in sorted(layer_dist.items()):
        print(f"  {layers} layers: {count} エントリ")

def main():
    """メイン実行関数"""
    print("🔍 Multi Base Pairing エントリ抽出スクリプト")
    print("=" * 60)
    
    # 入力・出力ファイル設定
    input_file = "analysis/pseudoknot_analysis_dssr_all.json"
    output_file = "analysis/multi_pairing_entries_dssr.json"
    
    # Step 1: エントリ抽出
    extracted_count = extract_multi_pairing_entries(input_file, output_file)
    
    # Step 2: 抽出したエントリを再読み込みして統計分析
    if extracted_count > 0:
        with open(output_file, 'r') as f:
            multi_pairing_entries = json.load(f)
        
        analyze_multi_pairing_stats(multi_pairing_entries)
    else:
        print("⚠️ Multi base pairingエントリが見つかりませんでした")

if __name__ == "__main__":
    main()
