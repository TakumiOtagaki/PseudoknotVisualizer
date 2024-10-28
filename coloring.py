from pymol import cmd
import json

def load_colors_from_json(file_path):
    """
    JSONファイルから深さごとの色設定を読み込む関数。
    
    :param file_path: JSONファイルのパス
    :return: 深さに対応する色の辞書
    """
    with open(file_path, 'r') as file:
        colors = json.load(file)
    return colors


def get_color_for_depth(depth, color_dict):
    """
    指定された深さに対応する色を取得する関数。
    
    :param depth: 深さ
    :param color_dict: 深さに対応する色の辞書
    :return: 指定された深さに対応する色（存在しない場合はデフォルト色を返す）
    """

    return color_dict.get(str(depth), color_dict.get("default"))


def coloring_canonical(pdb_object, chain, resi_i, color):
    # pdb_object の chain に対して、resi_i と resi_j の塩基を color で色付けする関数。
    cmd.color(color, f"{pdb_object} and chain {chain} and resi {resi_i}")
    return 
