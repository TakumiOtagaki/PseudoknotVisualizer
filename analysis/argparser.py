"""
Argument Parser for Pseudoknot Layer Analysis

Author: PseudoknotVisualizer Analysis
Date: 2025年7月21日
"""

import argparse


def create_parser():
    """
    コマンドライン引数パーサーを作成
    
    Returns:
        argparse.ArgumentParser: 設定済みのパーサー
    """
    parser = argparse.ArgumentParser(
        description="Pseudoknot Layer Analysis for RNA structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # DSSRを使用して解析
  python analysis/pseudoknotlayer_analysis.py --parser DSSR
  
  # RNAViewを使用して解析
  python analysis/pseudoknotlayer_analysis.py --parser RNAView
        """
    )
    
    # パーサー選択オプション
    parser.add_argument(
        "--parser",
        "-p",
        choices=["RNAView", "DSSR"],
        default="DSSR",
        help="RNA structure parser to use (default: DSSR)"
    )
    parser.add_argument(
        "--canonical-only",
        "-c",
        action="store_true",
        default=True,
        help="Only analyze canonical base pairs (default: False)"
    )
    
    return parser


def parse_args():
    """
    コマンドライン引数を解析
    
    Returns:
        argparse.Namespace: 解析された引数
    """
    parser = create_parser()
    return parser.parse_args()