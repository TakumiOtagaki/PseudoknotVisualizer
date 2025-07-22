"""
I/O Utilities for Pseudoknot Analysis

PDBファイルの列挙、チェーンID抽出、RNAView/DSSR実行を担当

Author: PseudoknotVisualizer Analysis
Date: 2025年7月21日
"""

import os
import sys
from pathlib import Path

# プロジェクトルートをパスに追加
script_dir = Path(__file__).parent.parent
sys.path.insert(0, str(script_dir))

from CLI_PseudoknotVisualizer import CLI_rnaview, CLI_dssr


def get_pdb_files(dataset_dir):
    """
    データセットディレクトリから全PDBファイルのリストを取得
    
    Args:
        dataset_dir (str): データセットディレクトリのパス
        
    Returns:
        list: PDBファイルのPathオブジェクトのリスト
    """
    pdb_files = list(Path(dataset_dir).glob("*.pdb"))
    print(f"Found {len(pdb_files)} PDB files in dataset")
    return sorted(pdb_files)


def extract_chain_from_filename(filename):
    """
    ファイル名からチェーン情報を抽出
    
    Args:
        filename (str): PDBファイル名
        
    Returns:
        str: チェーンID
        
    Examples:
        1EHZ_1_A.pdb -> A
        1J7T_1_A-B.pdb -> A (最初のチェーンを使用)
    """
    stem = Path(filename).stem
    parts = stem.split('_')
    if len(parts) >= 3:
        chain_part = parts[2]
        # ハイフンで区切られている場合は最初のチェーンを使用
        chain = chain_part.split('-')[0]
        return chain
    return 'A'  # デフォルト


def extract_actual_chain_from_pdb(pdb_file_path):
    """
    PDBファイルのREMARK 350行から実際のチェーンIDを抽出
    
    Args:
        pdb_file_path (str or Path): PDBファイルのパス
        
    Returns:
        str: チェーンID（見つからない場合は'A'）
    """
    try:
        with open(pdb_file_path, 'r') as f:
            for line in f:
                # REMARK 350 APPLY THE FOLLOWING TO CHAINS: の行を探す
                if line.startswith('REMARK 350') and 'APPLY THE FOLLOWING TO CHAINS:' in line:
                    # "CHAINS:" の後の部分を抽出
                    chains_part = line.split('CHAINS:')[1].strip()
                    # 最初のチェーンIDを取得（複数ある場合はスペースやカンマで区切られている）
                    chain_id = chains_part.split()[0].split(',')[0].strip()
                    return chain_id
                # ATOMレコードが始まったらヘッダー部分は終了
                elif line.startswith('ATOM'):
                    break
        
        # 見つからない場合はデフォルト
        print(f"Warning: No REMARK 350 CHAINS found in {pdb_file_path}, using default 'A'")
        return 'A'
        
    except Exception as e:
        print(f"Error reading PDB file {pdb_file_path}: {e}")
        return 'A'


def run_rnaview_analysis(pdb_file_path, chain_id):
    """
    RNAViewを実行して出力ファイルのパスを返す
    
    Args:
        pdb_file_path (str or Path): PDBファイルのパス
        chain_id (str): チェーンID
        
    Returns:
        Path: RNAView出力ファイルのパス（存在しない場合はNone）
    """
    print("HELLOOO")
    pdb_file = Path(pdb_file_path)
    print(f"Running RNAView for {pdb_file.name} with chain {chain_id}...")
    
    # RNAViewを実行
    raw_df = CLI_rnaview(str(pdb_file), chain_id)
    print(f"RNAView output generated:\n {raw_df}")
    # 出力ファイルパスを構築
    output_file = Path(f"intermediate/{pdb_file.name}.out")
    if not output_file.exists():
        raise FileNotFoundError(f"RNAView output not found for {pdb_file.name}")
    print(f"RNAView output generated: \n{output_file}")
    return output_file, raw_df


def run_dssr_analysis(pdb_file_path, chain_id):
    """
    DSSRを実行して出力ファイルのパスを返す
    
    Args:
        pdb_file_path (str or Path): PDBファイルのパス
        chain_id (str): チェーンID
        
    Returns:
        Path: DSSR出力ファイルのパス（存在しない場合はNone）, raw_df: pd.DataFrame
    """
    pdb_file = Path(pdb_file_path)
    print(f"Running DSSR for {pdb_file.name} with chain {chain_id}...")
    
    # DSSRを実行
    raw_df = CLI_dssr(str(pdb_file), chain_id)
    print(f"DSSR output generated: {raw_df}")
    
    # 出力ファイルパスを構築
    output_file = Path(f"intermediate/{pdb_file.name}.dssr.json")
    
    if not output_file.exists():
        print(f"Warning: DSSR output not found for {pdb_file}")
        raise FileNotFoundError(f"DSSR output not found for {pdb_file.name}")
    print(f"DSSR output generated: \n{output_file}")
    print(f"Warning: DSSR output not found for {pdb_file.name}")
    return output_file, raw_df


def run_parser_analysis(pdb_file_path, chain_id, parser="RNAView"):
    """
    指定されたパーサーを実行して出力ファイルのパスを返す
    
    Args:
        pdb_file_path (str or Path): PDBファイルのパス
        chain_id (str): チェーンID
        parser (str): "RNAView" or "DSSR"
        
    Returns:
        Path: 出力ファイルのパス（存在しない場合はNone）
    """
    if parser.upper() == "RNAVIEW":
        output_file, raw_df = run_rnaview_analysis(pdb_file_path, chain_id)
    elif parser.upper() == "DSSR":
        output_file, raw_df = run_dssr_analysis(pdb_file_path, chain_id)
    return output_file, raw_df
