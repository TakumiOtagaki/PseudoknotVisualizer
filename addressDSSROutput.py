import json
import pandas as pd
import re

def load_dssr_data(input_file: str):
    """
    DSSR JSON出力からベースペア情報をDataFrameとして読み込む
    
    Args:
        input_file (str): DSSRのJSON出力ファイルパス
        
    Returns:
        pd.DataFrame: 全てのベースペア情報
                     カラム: ["nt1", "nt2", "chain1", "chain2", "left_resi", "left_idx", 
                            "right_resi", "right_idx", "bp", "name", "saenger", "LW", "DSSR"]
    """
    with open(input_file, 'r') as infile:
        data = json.load(infile)

    result_data = []
    
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
            left_resi = nt1_match.group(2)  # 残基名 (G, C, A, U)
            left_idx = int(nt1_match.group(3))  # 残基番号
            right_resi = nt2_match.group(2)  # 残基名
            right_idx = int(nt2_match.group(3))  # 残基番号
            
            result_data.append({
                "nt1": pair["nt1"],
                "nt2": pair["nt2"], 
                "chain1": chain1,
                "chain2": chain2,
                "left_resi": left_resi,
                "left_idx": left_idx,
                "right_resi": right_resi,
                "right_idx": right_idx,
                "bp": pair.get("bp", ""),
                "name": pair.get("name", ""),
                "saenger": pair.get("Saenger", ""),
                "LW": pair.get("LW", ""),
                "DSSR": pair.get("DSSR", "")
            })

    # DataFrame に変換
    df = pd.DataFrame(result_data)
    return df


def canonical_extraction_from_dssr_df(df: pd.DataFrame):
    """
    DSSR DataFrameからカノニカルベースペアのみを抽出してフィルタリング
    
    Args:
        df (pd.DataFrame): load_dssr_data()で取得したDataFrame
        
    Returns:
        pd.DataFrame: フィルタリングされたベースペア情報
                     カラム: ["left_resi", "left_idx", "right_resi", "right_idx"]
    """
    if df.empty:
        return pd.DataFrame(columns=["left_resi", "left_idx", "right_resi", "right_idx"])
    
    # 同一チェーン内のベースペアのみを考慮
    same_chain_df = df[df["chain1"] == df["chain2"]]
    
    # Saenger番号によるカノニカルベースペアの判定
    canonical_df = same_chain_df[same_chain_df["saenger"].isin(["19-XIX", "20-XX", "28-XXVIII"])]
    
    # 必要なカラムのみを選択
    result_df = canonical_df[["left_resi", "left_idx", "right_resi", "right_idx"]].copy()
    
    return result_df


def extract_canonicalbp_from_dssr(df: pd.DataFrame):
    """
    DSSR JSON出力からカノニカルベースペア情報を抽出し、フィルタリングする
    （従来のextract_base_pairs_from_dssr関数の置き換え）
    
    Args:
        input_file (str): DSSRのJSON出力ファイルパス
        
    Returns:
        pd.DataFrame: フィルタリングされたベースペア情報
                     カラム: ["left_resi", "left_idx", "right_resi", "right_idx"]
    """
    # データをロード
    # df = load_dssr_data(input_file)
    # print(f"data loaded; Total base pairs found: {len(df)}")
    
    # カノニカルベースペアを抽出
    canonical_df = canonical_extraction_from_dssr_df(df)
    # print(f"canonical base pairs extracted; Total found: {len(canonical_df)}")
    
    return canonical_df


# 後方互換性のため、元の関数名も残す
# def extract_base_pairs_from_dssr(input_file: str):
#     """後方互換性のための関数（extract_canonicalbp_from_dssrの別名）"""
#     return extract_canonicalbp_from_dssr(input_file)


if __name__ == "__main__":
    dssr_output = "test/1KPD.dssr.json"
    
    # データローディングのテスト
    print("=== Testing data loading ===")
    full_df = load_dssr_data(dssr_output)
    print(f"Total base pairs found: {len(full_df)}")
    print(full_df.head())
    
    print("\n=== Testing canonical extraction ===")
    canonical_df = canonical_extraction_from_dssr_df(full_df)
    print(f"Canonical base pairs found: {len(canonical_df)}")
    print(canonical_df)
    
    print("\n=== Testing integrated function ===")
    integrated_df = extract_canonicalbp_from_dssr(dssr_output)
    print(f"Integrated function result: {len(integrated_df)}")
    print(integrated_df)
