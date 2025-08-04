from pymol import cmd

def viewer(pdb_object, chain_id, color="red"):
    """
    指定したオブジェクト内の特定のチェーンを赤色にする関数。
    
    :param pdb_object: PyMOL内のオブジェクト名
    :param chain_id: 対象のチェーンID
    :param color: 設定する色 (デフォルトは "red")
    """
    # オブジェクト内のチェーンを取得
    chains = cmd.get_chains(pdb_object)
    
    # 指定したチェーンが存在するか確認
    if chain_id not in chains:
        print(f"オブジェクト '{pdb_object}' にチェーン '{chain_id}' が存在しません。")
        return

    # 指定したチェーンを選択して色付け
    selection = f"{pdb_object} and chain {chain_id}"
    cmd.color(color, selection)
    print(f"{pdb_object} のチェーン {chain_id} を {color} に設定しました。")

# PyMOLのコマンドとして登録
cmd.extend("viewer", viewer)