"""
Parsers for RNA Structure Analysis

RNAView/DSSRの出力データを読み込み、共通フォーマットに変換

Author: PseudoknotVisualizer Analysis
Date: 2025年7月21日
"""

import sys
from pathlib import Path
import pandas as pd

# プロジェクトルートをパスに追加
script_dir = Path(__file__).parent.parent
sys.path.insert(0, str(script_dir))

from addressRNAviewOutput import load_rnaview_data, canonical_extraction_from_rnaview_df
from addressDSSROutput import load_dssr_data, raw_df_processing



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

def df_processing(df: pd.Dataframe, parser_type: str):
    """
    DataFrameから詳細な塩基対情報を共通フォーマットで作成
    
    Args:
        df (pandas.DataFrame): 塩基対データのDataFrame
        parser_type (str): "RNAView" or "DSSR"
        
    Returns:
        list: 塩基対詳細情報のリスト
    """
    bp_details = []
    for _, row in df.iterrows():
        if parser_type.upper() == "RNAVIEW":
            # RNAViewのSaenger分類を判定
            saenger = row.get("saenger", "")
            is_canonical = saenger in ["XX", "XIX", "XXVIII"]
            # wc_type = "Watson-Crick" if saenger in ["XIX", "XX"] else "Wobble" if saenger == "XXVIII" else "Non-WC"
            
        elif parser_type.upper() == "DSSR":
            # DSSRのSaenger分類を判定
            saenger = row.get("saenger", "")
            is_canonical = saenger in ["19-XIX", "20-XX", "28-XXVIII"]
            # wc_type = "Watson-Crick" if saenger in ["19-XIX", "20-XX"] else "Wobble" if saenger == "28-XXVIII" else "Non-WC"
        else:
            is_canonical = False
            # wc_type = "Unknown"
        
        bp_details.append({
            "position": [int(row["left_idx"]), int(row["right_idx"])],
            "residues": [row["left_resi"], row["right_resi"]],
            "is_canonical": is_canonical,
            # "wc_type": wc_type,
            "saenger_id": saenger
        })
    return bp_details


def parse_rnaview_output(output_file_path):
    """
    RNAView出力ファイルを解析して共通フォーマットで返す
    
    Args:
        output_file_path (str or Path): RNAView出力ファイルのパス
        
    Returns:
        tuple: (全塩基対詳細情報, canonical塩基対詳細情報, canonical塩基対リスト)
    """
    # データをロード
    all_bp_df = load_rnaview_data(str(output_file_path))
    # canonical_bp_df = canonical_extraction_from_rnaview_df(all_bp_df)
    all_bp_df = raw_df_processing(all_bp_df)
    
    # 詳細情報を作成
    all_bp_details = create_bp_details(all_bp_df, "RNAView")
    canonical_bp_details = create_bp_details(canonical_bp_df, "RNAView")
    
    # PKextractor用の位置情報のみのリスト
    canonical_bp_list = [(row["left_idx"], row["right_idx"]) for _, row in canonical_bp_df.iterrows()]
    
    return all_bp_details, canonical_bp_details, canonical_bp_list




def parse_output_file(output_file_path, parser_type):
    """
    指定されたパーサーの出力ファイルを解析して共通フォーマットで返す
    
    Args:
        output_file_path (str or Path): 出力ファイルのパス
        parser_type (str): "RNAView" or "DSSR"
        
    Returns:
        tuple: (全塩基対詳細情報, canonical塩基対詳細情報, canonical塩基対リスト)
    """
    if parser_type.upper() == "RNAVIEW":
        raw_df = parse_rnaview_output(output_file_path)
    elif parser_type.upper() == "DSSR":
        raw_df = load_dssr_data(str(output_file_path))
    processed_df = raw_df_processing(raw_df, "DSSR")
    return processed_df


def filter_abnormal_pairs(bp_details):
    """
    自己ペア（i=j）を除外したリストを作成
    
    Args:
        bp_details (list): 塩基対詳細情報のリスト
        bp_list (list): 塩基対位置のリスト
        
        - (i, i) のような self pairing を除外
        - (i, j) と (i, j') のような塩基対があれば、一方が canonical base pair ならば、そちらを残す。
            - 両方とも non-canonical ならば両方とも無視する。
            - 両方とも canonical ならば、raise Error
    Returns:
        tuple: (basepair details without abnormal pairs, abnormal pairs list)

    """
    # 自己ペア（i=j）を検出
    abnormal_pairs = [bp["position"] for bp in bp_details if bp["position"][0] == bp["position"][1]]
    bp_details_filtered = [bp for bp in bp_details if bp["position"][0] != bp["position"][1]]


    # (i, j) と (i, j') のような塩基対を検出
    for i, bp1 in enumerate(bp_details):
        for j in range(i + 1, len(bp_details)):
            bp2 = bp_details[j]
            if bp1["position"][0] == bp2["position"][0] and bp1["position"][1] == bp2["position"][1]:
                # 同じ位置のペアが見つかった
                if bp1["is_canonical"] and not bp2["is_canonical"]:
                    abnormal_pairs.append(bp1)
                elif not bp1["is_canonical"] and bp2["is_canonical"]:
                    abnormal_pairs.append(bp2)
                elif bp1["is_canonical"] and bp2["is_canonical"]:
                    raise ValueError("Both pairs are canonical")
                else:
                    # 両方とも non-canonical ならば無視
                    continue

    return bp_details_filtered, self_pairs, abnormal_pairs
