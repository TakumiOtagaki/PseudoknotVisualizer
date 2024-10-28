import re
import pandas as pd

# Only canonical base pairs are considered!!
# Considering non-canonical base pairs is a future work. 2024.10.28
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
        # print(f"fields: {fields}")
        
        # 特定条件のチェック
        if fields[1] == fields[5]:  # considering only chain-internal bp
            if len(fields) != 9: continue
            if fields[8] in ["XX", "XIX", "XXVIII"]:  # Canonical BP
                # result_lines.append("\t".join(fields[:9]))
                # print("hi")

                left_idx, right_idx = tuple(map(int, fields[0].replace(",", "").split("_")))
                left_resi, right_resi = tuple(map(str, fields[3].split("-")))
                
                # print(f"left_resi={left_resi}, right_resi={right_resi}, left_idx={left_idx}, right_idx={right_idx}")
                result_lines.append((left_resi, left_idx, right_resi, right_idx))

                
    # フィルタリング結果を pandas DataFrame に変換
    df = pd.DataFrame(result_lines, columns=["left_resi", "left_idx", "right_resi", "right_idx"])
    return df


if __name__ == "__main__":
    rnaview_output = "/large/otgk/PseudoknotVisualizer/intermediate/1KPD.pdb.out"
    df = extract_base_pairs(rnaview_output)
    print(df)

