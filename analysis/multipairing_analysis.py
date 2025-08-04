#!/usr/bin/env python3
"""
Multi-pairing Analysis Script

JSONファイル生成時にmulti-pairing（triplet等）がどう処理されているかを確認
"""

import json
from pathlib import Path

def analyze_multipairing_handling():
    """Multi-pairing処理の確認"""
    
    # サンプルファイルを調査
    json_file = "analysis/pseudoknot_analysis_dssr_all.json"
    
    if not Path(json_file).exists():
        print(f"❌ {json_file} が見つかりません")
        return
    
    with open(json_file, 'r') as f:
        results = json.load(f)
    
    print("=" * 60)
    print("Multi-pairing (Triplet) 処理分析")
    print("=" * 60)
    
    # 重複canonical塩基対が多いチェーンを探す
    high_dup_chains = []
    for result in results:
        dup_pairs = result.get('dup_canonical_pairs', [])
        if len(dup_pairs) > 10:  # 10個以上の重複があるチェーン
            high_dup_chains.append({
                'pdb_id': result['pdb_id'],
                'chain_id': result['chain_id'],
                'dup_count': len(dup_pairs),
                'total_bp': result.get('total_bp_count', 0),
                'abnormal_count': len(result.get('abnormal_pairs', [])),
                'dup_pairs': dup_pairs[:5]  # 最初の5個だけ表示
            })
    
    print(f"重複canonical塩基対が10個以上あるチェーン: {len(high_dup_chains)}個")
    print()
    
    # 上位5つを詳細表示
    high_dup_chains.sort(key=lambda x: x['dup_count'], reverse=True)
    for i, chain in enumerate(high_dup_chains[:5]):
        print(f"🔍 #{i+1}: {chain['pdb_id']} (Chain {chain['chain_id']})")
        print(f"   重複canonical: {chain['dup_count']}個")
        print(f"   総塩基対: {chain['total_bp']}個")
        print(f"   異常ペア: {chain['abnormal_count']}個")
        print(f"   重複例: {chain['dup_pairs']}")
        print()
    
    # 処理方法の確認
    print("📋 現在の処理方法:")
    print("1. Multi-pairing検出時:")
    print("   - 同じ塩基位置に複数の塩基対がある場合を検出")
    print("   - canonical vs non-canonical → canonicalを保持")
    print("   - canonical vs canonical → 両方をdup_canonical_pairsに記録")
    print("   - non-canonical vs non-canonical → 両方を除外")
    print()
    print("2. 最終的なJSONファイルには:")
    print("   - フィルタリング後の塩基対のみ記録")
    print("   - 除外された塩基対は abnormal_pairs に記録")
    print("   - 重複canonical塩基対は dup_canonical_pairs に記録")
    print("   - エントリ全体は除外せず、問題のある塩基対のみ除外")

if __name__ == "__main__":
    analyze_multipairing_handling()
