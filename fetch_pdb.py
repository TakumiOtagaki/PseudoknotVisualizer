import urllib.request

def download_pdb(pdb_id, filename=None):
    """
    指定された PDB ID を参照して、PDB ファイルをダウンロードし、
    given_ID.pdb として保存します。
    
    Parameters:
    pdb_id (str): PDB の ID（例: '1ABC'）
    filename (str): 保存するファイル名（指定しない場合は given_ID.pdb として保存）
    """
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    if not filename:
        filename = f"{pdb_id}.pdb"
    
    try:
        urllib.request.urlretrieve(url, filename)
        print(f"PDB file for {pdb_id} downloaded as {filename}")
    except Exception as e:
        print(f"Failed to download PDB file for {pdb_id}: {e}")

# 使用例
pdb_id = input("Enter PDB ID: ")
download_pdb(pdb_id, filename="given_ID.pdb")