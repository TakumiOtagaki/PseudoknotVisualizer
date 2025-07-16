import json
import pandas as pd
import re

# Only canonical base pairs are considered!!
# Considering non-canonical base pairs is a future work. 2024.10.28
def extract_base_pairs_from_dssr(input_file: str):
    """
    DSSR JSON出力からベースペア情報を抽出し、フィルタリングする
    
    Args:
        input_file (str): DSSRのJSON出力ファイルパス
        
    Returns:
        pd.DataFrame: フィルタリングされたベースペア情報
                     カラム: ["left_resi", "left_idx", "right_resi", "right_idx"]
    """
    with open(input_file, 'r') as infile:
        data = json.load(infile)

    result_lines = []
    
    # pairsセクションが存在する場合のみ処理
    if "pairs" in data:
        for pair in data["pairs"]:
            # nt1とnt2から チェーン、残基名、残基番号を抽出
            # 例: "A.G1" -> チェーン="A", 残基="G", 番号="1"
            nt1_match = re.match(r"([A-Z]+)\.([A-Z])(\d+)", pair["nt1"])
            nt2_match = re.match(r"([A-Z]+)\.([A-Z])(\d+)", pair["nt2"])
            
            if not nt1_match or not nt2_match:
                continue
                
            # チェーン情報を取得
            chain1 = nt1_match.group(1)
            chain2 = nt2_match.group(1)
            
            # 同一チェーン内のベースペアのみを考慮
            if chain1 == chain2:
                # Saenger番号によるカノニカルベースペアの判定
                saenger = pair.get("Saenger", "")
                if saenger in ["19-XIX", "20-XX", "28-XXVIII"]:  # Canonical BP
                    left_resi = nt1_match.group(2)  # 残基名 (G, C, A, U)
                    left_idx = int(nt1_match.group(3))  # 残基番号
                    right_resi = nt2_match.group(2)  # 残基名
                    right_idx = int(nt2_match.group(3))  # 残基番号
                    
                    result_lines.append((left_resi, left_idx, right_resi, right_idx))

    # フィルタリング結果を pandas DataFrame に変換
    df = pd.DataFrame(result_lines, columns=["left_resi", "left_idx", "right_resi", "right_idx"])
    return df


if __name__ == "__main__":
    dssr_output = "/Users/ootagakitakumi/PseudoknotVisualizer/test/1KPD.dssr.json"
    df = extract_base_pairs_from_dssr(dssr_output)
    print(df)
    print(f"Number of canonical base pairs found: {len(df)}")
