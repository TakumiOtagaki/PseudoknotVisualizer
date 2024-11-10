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
pdb_id = input("Enter PDB ID (q to quit): ")
if pdb_id == "q":
    exit()
output_file = input("Enter output filename(if not provided, pdb_id.pdb will be created in current directory): ")
if output_file == "":
    output_file = "./" + pdb_id + ".pdb"
download_pdb(pdb_id, output_file)