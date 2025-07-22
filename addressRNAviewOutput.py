import re
import pandas as pd

def load_rnaview_data(input_file: str):
    """
    RNAView出力からベースペア情報をDataFrameとして読み込む
    
    Args:
        input_file (str): RNAViewの出力ファイルパス
        
    Returns:
        pd.DataFrame: 全てのベースペア情報
                     カラム: ["pair_idx", "chain1", "chain2", "left_resi", "left_idx", 
                            "right_resi", "right_idx", "saenger", "additional_fields"]
    """
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # 'BEGIN_base-pair' と 'END_base-pair' の間のテーブル行を抽出
    extracting = False
    table_lines = []
    for line in lines:
        if "BEGIN_base-pair" in line:
            extracting = True
            continue
        if "END_base-pair" in line:
            extracting = False
            break
        if extracting:
            table_lines.append(line.strip())

    result_data = []
    for line in table_lines:
        # タブで区切られたフィールドに変換
        fields = re.split(r"\s+", line)
        
        if len(fields) < 6:  # 最低限必要なフィールド数をチェック
            continue
            
        try:
            # インデックス情報の抽出
            left_idx, right_idx = tuple(map(int, fields[0].replace(",", "").split("_")))
            chain1 = fields[1]
            chain2 = fields[5]
            
            # 残基名の抽出
            left_resi, right_resi = tuple(map(str, fields[3].split("-")))
            
            # Saenger番号（最後のフィールド）
            saenger = fields[-1] if len(fields) > 0 else ""
            
            result_data.append({
                "pair_idx": fields[0],
                "chain1": chain1,
                "chain2": chain2,
                "left_resi": left_resi,
                "left_idx": left_idx,
                "right_resi": right_resi,
                "right_idx": right_idx,
                "saenger": saenger,
                "additional_fields": fields[2:],  # その他のフィールド
            })
        except (ValueError, IndexError) as e:
            # パースに失敗した行はスキップ
            continue

    # DataFrame に変換
    df = pd.DataFrame(result_data)
    return df


def canonical_extraction_from_rnaview_df(df: pd.DataFrame):
    """
    RNAView DataFrameからカノニカルベースペアのみを抽出してフィルタリング
    
    Args:
        df (pd.DataFrame): load_rnaview_data()で取得したDataFrame
        
    Returns:
        pd.DataFrame: フィルタリングされたベースペア情報
                     カラム: ["left_resi", "left_idx", "right_resi", "right_idx"]
    """
    if df.empty:
        return pd.DataFrame(columns=["left_resi", "left_idx", "right_resi", "right_idx"])
    
    # 同一チェーン内のベースペアのみを考慮
    same_chain_df = df[df["chain1"] == df["chain2"]]
    
    # Saenger番号によるカノニカルベースペアの判定
    canonical_df = same_chain_df[same_chain_df["saenger"].isin(["XX", "XIX", "XXVIII"])]
    
    # 必要なカラムのみを選択
    result_df = canonical_df[["left_resi", "left_idx", "right_resi", "right_idx"]].copy()
    
    return result_df


def extract_canonicalbp_from_rnaview(df: pd.DataFrame):
    """
    RNAView出力からカノニカルベースペア情報を抽出し、フィルタリングする
    （従来のextract_base_pairs_from_rnaview関数の置き換え）
    
    Args:
        input_file (str): RNAViewの出力ファイルパス
        
    Returns:
        pd.DataFrame: フィルタリングされたベースペア情報
                     カラム: ["left_resi", "left_idx", "right_resi", "right_idx"]
    """
    # データをロード
    # df = load_rnaview_data(input_file)
    
    # カノニカルベースペアを抽出
    canonical_df = canonical_extraction_from_rnaview_df(df)
    
    return canonical_df


# 後方互換性のため、元の関数名も残す
def extract_base_pairs_from_rnaview(input_file: str):
    """後方互換性のための関数（extract_canonicalbp_from_rnaviewの別名）"""
    return extract_canonicalbp_from_rnaview(input_file)


if __name__ == "__main__":
    rnaview_output = "test/1KPD.pdb.out"
    
    # データローディングのテスト
    print("=== Testing data loading ===")
    full_df = load_rnaview_data(rnaview_output)
    print(f"Total base pairs found: {len(full_df)}")
    print(full_df.head())
    
    print("\n=== Testing canonical extraction ===")
    canonical_df = canonical_extraction_from_rnaview_df(full_df)
    print(f"Canonical base pairs found: {len(canonical_df)}")
    print(canonical_df)
    
    print("\n=== Testing integrated function ===")
    integrated_df = extract_canonicalbp_from_rnaview(rnaview_output)
    print(f"Integrated function result: {len(integrated_df)}")
    print(integrated_df)
