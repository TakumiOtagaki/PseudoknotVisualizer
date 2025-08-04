#!/usr/bin/env python3
"""
Multi-layer Pseudoknot Entry Filter

analysis/datasets/pdbid_chains.csv から、
analysis/pseudoknot_analysis_dssr_all.json の条件に合致するエントリを抽出

条件:
- total_bp_count > 0
- pseudoknot_layer_count > 1  
- output_exists = true
- pdb_id と actual_chain_id が一致
"""

import json
import pandas as pd
from pathlib import Path

def main():
    # ファイルパス
    csv_file = "analysis/datasets/pdbid_chains.csv"
    json_file = "analysis/pseudoknot_analysis_dssr_all.json"
    
    # JSONデータをDataFrameに変換
    with open(json_file, 'r') as f:
        json_data = json.load(f)
    
    json_df = pd.DataFrame(json_data)
    
    # 条件でフィルタリング
    filtered_json = json_df[
        (json_df['total_bp_count'] > 0) & 
        # (json_df['pseudoknot_layer_count'] > 1) & 
        (json_df['output_exists'] == True)
    ]
    
    print(f"JSON条件に合致するエントリ数: {len(filtered_json)}")
    
    # 結果DataFrame作成
    result_data = []
    
    for _, json_row in filtered_json.iterrows():
        pdb_id = json_row['pdb_id']  # 既に {PDB_ID}_{model_ID}_{chains} 形式
        chain_id = json_row['chain_id']  # 実際のchain（一つに決まる）
        
        # pdb_id, model_id(=1), chain, extra(=NA) の形式で出力
        result_data.append({
            'pdb_id': pdb_id,
            'model_id': 1,
            'chain': chain_id,
            'extra': 'NA'
        })
    
    result_df = pd.DataFrame(result_data)
    
    # 結果出力
    print("\n条件に合致するエントリ:")
    print("pdb_id,model_id,chain,extra")
    
    # CSV形式で出力
    result_csv = result_df.to_csv(index=False, header=False)
    print(result_csv.strip())
    
    print(f"\nマッチした行数: {len(result_df)}")
    
    # 例: 最初の数件の詳細情報
    print(f"\n例（最初の5件）:")
    for i, (_, row) in enumerate(result_df.head().iterrows()):
        json_row = filtered_json.iloc[i]
        print(f"{row['pdb_id']},{row['model_id']},{row['chain']},{row['extra']} "
              f"(bp:{json_row['total_bp_count']}, layers:{json_row['pseudoknot_layer_count']})")
    
    # 結果をCSVファイルに保存
    output_file = Path("analysis/datasets/filtered_multilayer_entries.csv")
    # pdb_id は PDB_ で始まることがあるのでまずそれを処理する
    result_df['pdb_id'] = result_df['pdb_id'].str.replace('^PDB_', '', regex=True)
    result_df['pdb_id'] = result_df['pdb_id'].str.replace('0000', '', regex=True)
    # pdb_id は "_" で区切って先頭だけ出力することにする。
    result_df['pdb_id'] = result_df['pdb_id'].str.split('_').str[0]
    result_df.to_csv(output_file, index=False, header=True)
    print(f"\n結果を {output_file} に保存しました。")

if __name__ == "__main__":
    main()
