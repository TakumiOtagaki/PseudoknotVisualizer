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
    python analysis/pseudoknotlayer_analysis.py --annotator DSSR
  
    # RNAViewを使用して解析
    python analysis/pseudoknotlayer_analysis.py --annotator RNAView
        """
    )
    
    # パーサー選択オプション
    parser.add_argument(
        "--annotator",
        "-a",
        choices=["RNAView", "DSSR"],
        default="DSSR",
        help="RNA structure annotator to use (default: DSSR)"
    )
    # Hidden legacy flags for backward compatibility
    parser.add_argument("--parser", dest="annotator", choices=["RNAView", "DSSR"], help=argparse.SUPPRESS)
    parser.add_argument("-p", dest="annotator", choices=["RNAView", "DSSR"], help=argparse.SUPPRESS)
    parser.add_argument(
        "--canonical-only",
        "-c",
        action="store_true",
        default=False,
        help="Only analyze canonical base pairs (default: False)"
    )
    parser.add_argument(
        "--ncpus",
        "-n",
        type=int,
        default=1,
        help="Number of CPU cores to use for parallel processing (default: 1)"
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