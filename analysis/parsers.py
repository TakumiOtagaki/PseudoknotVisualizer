"""
Parsers for RNA Structure Analysis

RNAView/DSSRの出力データを読み込み、共通フォーマットに変換

Author: PseudoknotVisualizer Analysis
Date: 2025年7月21日
"""

import sys
from pathlib import Path

# プロジェクトルートをパスに追加
script_dir = Path(__file__).parent.parent
sys.path.insert(0, str(script_dir))

from addressRNAviewOutput import load_rnaview_data, canonical_extraction_from_rnaview_df
from addressDSSROutput import load_dssr_data, canonical_extraction_from_dssr_df


def create_bp_details(df, parser_type):
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
            wc_type = "Watson-Crick" if saenger in ["XIX", "XX"] else "Wobble" if saenger == "XXVIII" else "Non-WC"
            
        elif parser_type.upper() == "DSSR":
            # DSSRのSaenger分類を判定
            saenger = row.get("saenger", "")
            is_canonical = saenger in ["19-XIX", "20-XX", "28-XXVIII"]
            wc_type = "Watson-Crick" if saenger in ["19-XIX", "20-XX"] else "Wobble" if saenger == "28-XXVIII" else "Non-WC"
        else:
            is_canonical = False
            wc_type = "Unknown"
        
        bp_details.append({
            "position": [int(row["left_idx"]), int(row["right_idx"])],
            "residues": [row["left_resi"], row["right_resi"]],
            "is_canonical": is_canonical,
            "wc_type": wc_type,
            "saenger_id": saenger,
            "base_pair_type": "canonical" if is_canonical else "non-canonical"
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
    canonical_bp_df = canonical_extraction_from_rnaview_df(all_bp_df)
    
    # 詳細情報を作成
    all_bp_details = create_bp_details(all_bp_df, "RNAView")
    canonical_bp_details = create_bp_details(canonical_bp_df, "RNAView")
    
    # PKextractor用の位置情報のみのリスト
    canonical_bp_list = [(row["left_idx"], row["right_idx"]) for _, row in canonical_bp_df.iterrows()]
    
    return all_bp_details, canonical_bp_details, canonical_bp_list


def parse_dssr_output(output_file_path):
    """
    DSSR出力ファイルを解析して共通フォーマットで返す
    
    Args:
        output_file_path (str or Path): DSSR出力ファイルのパス
        
    Returns:
        tuple: (全塩基対詳細情報, canonical塩基対詳細情報, canonical塩基対リスト)
    """
    # データをロード
    all_bp_df = load_dssr_data(str(output_file_path))
    canonical_bp_df = canonical_extraction_from_dssr_df(all_bp_df)
    
    # 詳細情報を作成
    all_bp_details = create_bp_details(all_bp_df, "DSSR")
    canonical_bp_details = create_bp_details(canonical_bp_df, "DSSR")
    
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
        return parse_rnaview_output(output_file_path)
    elif parser_type.upper() == "DSSR":
        return parse_dssr_output(output_file_path)
    else:
        raise ValueError(f"Unsupported parser: {parser_type}")


def filter_self_pairs(bp_details, bp_list):
    """
    自己ペア（i=j）を除外したリストを作成
    
    Args:
        bp_details (list): 塩基対詳細情報のリスト
        bp_list (list): 塩基対位置のリスト
        
    Returns:
        tuple: (フィルタリング後の詳細情報, フィルタリング後の位置リスト, 自己ペアリスト)
    """
    # 自己ペア（i=j）を検出
    self_pairs = [bp for bp in bp_list if bp[0] == bp[1]]
    
    # 自己ペアを除外したリストを作成
    bp_details_filtered = [bp for bp in bp_details if bp["position"][0] != bp["position"][1]]
    bp_list_filtered = [bp for bp in bp_list if bp[0] != bp[1]]
    
    return bp_details_filtered, bp_list_filtered, self_pairs
