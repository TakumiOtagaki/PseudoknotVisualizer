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
    if not input_file.exists():
        return pd.DataFrame(columns=["nt1", "nt2", "chain1", "chain2", "left_resi", "left_idx", 
                                     "right_resi", "right_idx", "bp", "name", "saenger", "LW", "DSSR"])
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


# def raw_DSSR_df_processing(df: pd.DataFrame):
#     """
#     DSSRの生データを処理して、必要なカラムのみを抽出し、カノニカルベースペアのフラグを追加
#     samle chain におけるものだけを取り出すことにしている
#     """

#     if df.empty:
#         return pd.DataFrame(columns=["left_resi", "left_idx", "right_resi", "right_idx", "chain1", "chain2", "is_canonical", "saenger_id"])
    
#     # 同一チェーン内のベースペアのみを考慮
#     same_chain_df = df[df["chain1"] == df["chain2"]]
    
#     # Saenger番号によるカノニカルベースペアの判定
#     # processed_df = same_chain_df[same_chain_df["saenger"].isin(["19-XIX", "20-XX", "28-XXVIII"])]
#     processed_df = same_chain_df.copy()
#     processed_df["is_canonical"] = same_chain_df["saenger"].isin(["19-XIX", "20-XX", "28-XXVIII"])
    
#     # 必要なカラムのみを選択
#     result_df = processed_df[["left_resi", "left_idx", "right_resi", "right_idx", "chain1", "chain2", "is_canonical", "saenger"]]
#     return result_df
