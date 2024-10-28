import re
import pandas as pd

def extract_base_pairs(input_file: str):
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

    # 抽出した行をフィルタリング
    result_lines = []
    for line in table_lines:
        # タブで区切られたフィールドに変換
        fields = re.split(r"\s+", line)
        
        # 特定条件のチェック
        if len(fields) >= 10 and fields[1] == fields[5]:  # chain-internal bp
            if fields[8] in ["XX", "XIX", "XXVIII"]:  # Canonical BP
                result_lines.append("\t".join(fields[:9]))  # 最初の9列まで出力
                
    # フィルタリング結果を pandas DataFrame に変換
    df = pd.read_csv(
        pd.compat.StringIO("\n".join(result_lines)),
        sep="\t",
        names=[
            "residue1", "chain1", "residue2", "chain2", "orientation", "type", "edge1", "edge2", "distance"
        ]
    )
