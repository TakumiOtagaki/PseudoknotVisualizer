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

from addressRNAviewOutput import load_rnaview_data
from addressDSSROutput import load_dssr_data


def raw_df_processing(df: pd.DataFrame, parser_type: str):
    """
    DataFrameから詳細な塩基対情報を共通フォーマットで作成
    is_canonicalフラグを追加し、Saenger分類に基づいて塩基対を分類
    
    Args:
        df (pandas.DataFrame): 塩基対データのDataFrame
        parser_type (str): "RNAView" or "DSSR"
        
    Returns:
        list: 塩基対詳細情報のリスト
    """
    # print(f"DataFrame shape: {df.shape}")
    # print(f"Columns: {df.columns.tolist()}")
    # print(f"First few rows:\n{df.head()}")

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
            saenger = ""
        bp_details.append({
            "position": [int(row["left_idx"]), int(row["right_idx"])],
            "residues": [row["left_resi"], row["right_resi"]],
            "is_canonical": is_canonical,
            # "wc_type": wc_type,
            "saenger_id": saenger
        })
    # position の left について昇順にソート
    bp_details.sort(key=lambda x: x["position"][0])
    # 空の場合は空のDataFrameを返す
    raw_df = pd.DataFrame(bp_details)
    if raw_df.empty:
        return pd.DataFrame(columns=["position", "residues", "is_canonical", "saenger_id"])
    return pd.DataFrame(bp_details)


def parse_output_file(output_file_path, parser_type):
    """
    指定されたパーサーの出力ファイルを解析して共通フォーマットで返す
    
    Args:
        output_file_path (str or Path): 出力ファイルのパス
        parser_type (str): "RNAView" or "DSSR"
        
    Returns:
        tuple: 全塩基対詳細情報
    """
    if parser_type.upper() == "RNAVIEW":
        raw_df = load_rnaview_data(str(output_file_path))
    elif parser_type.upper() == "DSSR":
        raw_df = load_dssr_data(str(output_file_path))
    # processed_df = raw_df_processing(raw_df, parser_type)
    return raw_df


def filter_abnormal_pairs(processed_df: pd.DataFrame):
    """
    自己ペア（i=j）を除外したリストを作成
    
    Args:
        processed_df (pd.DataFrame): 塩基対詳細情報のDataFrame
        - (i, i) のような self pairing を除外
        - (i, j) と (i, j') のような塩基対があれば、一方が canonical base pair ならば、そちらを残す。
            - 両方とも non-canonical ならば両方とも無視する。
            - 両方とも canonical ならば、raise Error
    Returns:
        tuple: (basepair details without abnormal pairs, abnormal pairs list)

    """
    processed_dict = processed_df.to_dict(orient='records')
    # 自己ペア（i=j）を検出
    abnormal_pairs = [bp["position"] for bp in processed_dict if bp["position"][0] == bp["position"][1]]
    processed_dict_filtered = [bp for bp in processed_dict if bp["position"][0] != bp["position"][1]]

    base_pairs = processed_df["position"].values.tolist()
    base_pairs = [tuple(bp) for bp in base_pairs]
    for bp in processed_df.values:
        # print("hello")
        left_idx_tmp, right_idx_tmp = tuple(bp[0])
        left_idx = min(left_idx_tmp, right_idx_tmp)
        right_idx = max(left_idx_tmp, right_idx_tmp)
        # print(f"left_idx: {left_idx}, right_idx: {right_idx}")
        if (right_idx, left_idx) in base_pairs:
            print(f"Switched pair found: {(left_idx, right_idx)} and {(right_idx, left_idx)}")
            abnormal_pairs.append([right_idx, left_idx])
        # bp[position][0] と bp[position][1] を入れ替えて残りはそのままにしたものが存在していたら一方だけにする
        # bp_pos_switched = (bp["position"][1], bp["position"][0])
        # if bp_pos_switched in processed_dict_filtered:
        #     processed_dict_filtered.remove(bp_pos_switched)

    # (i, j) と (i, j') のような塩基対を検出
    dup_canonical_pairs = set()
    for i, bp1 in enumerate(processed_dict_filtered):
        for j in range(i + 1, len(processed_dict_filtered)):
            bp2 = processed_dict_filtered[j]
            if (set(bp1["position"]) & set(bp2["position"])):
                # 同じ位置のペアが見つかった
                if bp1["is_canonical"] and not bp2["is_canonical"]:
                    abnormal_pairs.append(bp2["position"])
                elif not bp1["is_canonical"] and bp2["is_canonical"]:
                    abnormal_pairs.append(bp1["position"])
                elif bp1["is_canonical"] and bp2["is_canonical"]:
                    # raise ValueError("Both pairs are canonical")
                    dup_canonical_pairs.add((bp1["position"][0], bp1["position"][1]))
                    dup_canonical_pairs.add((bp2["position"][0], bp2["position"][1]))
                else:
                    # 両方とも non-canonical ならば無視
                    abnormal_pairs.append(bp1["position"])
                    abnormal_pairs.append(bp2["position"])

    bp_details_filtered = [bp for bp in processed_dict_filtered if bp["position"] not in abnormal_pairs]
    bp_details_filtered_df = pd.DataFrame(bp_details_filtered)
    if bp_details_filtered_df.empty:
        return pd.DataFrame(columns=["position", "residues", "is_canonical", "saenger_id"]), abnormal_pairs, dup_canonical_pairs
    # print(f"Filtered out {len(abnormal_pairs)} abnormal pairs")
    print(f"Remaining pairs: {len(bp_details_filtered)}")
    return bp_details_filtered_df, abnormal_pairs, dup_canonical_pairs
