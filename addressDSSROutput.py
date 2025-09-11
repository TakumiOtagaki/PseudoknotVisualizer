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
    # Accept Path-like by converting to str where needed
    try:
        from pathlib import Path as _Path
        _p = _Path(input_file)
    except Exception:
        _p = None
    if _p is not None and not _p.exists():
        return pd.DataFrame(columns=["nt1", "nt2", "chain1", "chain2", "left_resi", "left_idx", 
                                     "right_resi", "right_idx", "bp", "name", "saenger", "LW", "DSSR"])
    with open(str(input_file), 'r') as infile:
        data = json.load(infile)

    result_data = []
    
    # pairsセクションが存在する場合のみ処理
    if "pairs" in data:
        for pair in data["pairs"]:
            # nt1/nt2 の形式は DSSR では一般に
            #   <chain>.<base><seqnum>[insertion]
            # のようになっている（例: "A.G1", 数字チェーンの例: "5.C1"、挿入コード付き: "A.G10A"）。
            # 既存の実装はチェーンIDを [A-Z]+ に縛っていたため、'5' などの数値チェーンで失敗していた。
            # より寛容なパターンでパースする。
            # - チェーンID: ドット以外の連続文字
            # - 塩基: 先頭1文字（A/C/G/U など）
            # - 残基番号: 任意桁の整数（負も許容）＋任意の挿入コード1文字（あれば捨てる）
            nt_pattern = r"([^\.]+)\.([A-Za-z])(-?\d+)([A-Za-z]?)"
            nt1_match = re.match(nt_pattern, pair.get("nt1", ""))
            nt2_match = re.match(nt_pattern, pair.get("nt2", ""))

            if not nt1_match or not nt2_match:
                # それでも合わないケースはスキップ（後続処理で空DFは安全に扱われる）
                continue

            # チェーン情報と残基情報
            chain1 = nt1_match.group(1)
            chain2 = nt2_match.group(1)
            left_resi = nt1_match.group(2)   # 残基名 (G, C, A, U など)
            right_resi = nt2_match.group(2)
            # 残基番号（整数部を使用。挿入コードは無視）
            left_idx = int(nt1_match.group(3))
            right_idx = int(nt2_match.group(3))

            result_data.append({
                "nt1": pair.get("nt1", ""),
                "nt2": pair.get("nt2", ""),
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
