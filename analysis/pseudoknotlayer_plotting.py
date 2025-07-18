"""
Pseudoknot Layer Plotting Script

このスクリプトはpseudoknotlayer_analysis.pyの結果を可視化します：
1. pseudoknot_analysis_rnaview.json と pseudoknot_analysis_dssr.json を読み込み
2. Main Layer (Layer 0) とPseudoknot Layers (Layer 1+) の統計を可視化
3. Canonical/Non-canonical base pairの分布を解析

Author: PseudoknotVisualizer Analysis
Date: 2025年7月18日
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# 日本語フォント設定
plt.rcParams['font.family'] = ['DejaVu Sans', 'Hiragino Sans', 'Yu Gothic', 'Meiryo', 'Takao', 'IPAexGothic', 'IPAPGothic', 'VL PGothic', 'Noto Sans CJK JP']

# スクリプトディレクトリ
script_dir = Path(__file__).parent
analysis_dir = script_dir

def load_analysis_results(parser_name):
    """
    解析結果JSONファイルを読み込み
    
    Args:
        parser_name (str): "rnaview" or "dssr"
        
    Returns:
        list: 解析結果のリスト
    """
    json_file = analysis_dir / f"pseudoknot_analysis_{parser_name.lower()}.json"
    
    if not json_file.exists():
        print(f"Warning: {json_file} not found!")
        return []
    
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    print(f"Loaded {len(data)} structures from {json_file}")
    return data

def extract_layer_statistics(results):
    """
    レイヤー統計情報を抽出
    
    Args:
        results (list): 解析結果のリスト
        
    Returns:
        dict: 統計情報
    """
    stats = {
        'main_layer_stats': [],      # Layer 0 (Main/Core Layer)
        'pseudoknot_layer_stats': [], # Layer 1+ (Pseudoknot Layers)
        'total_stats': [],           # 全体統計
        'structure_info': []         # 構造情報
    }
    
    for result in results:
        pdb_id = result['pdb_id']
        parser = result['parser']
        total_bp = result['total_bp_count']
        total_canonical = result['total_canonical_bp_count']
        total_layers = result['pseudoknot_layer_count']
        
        # 全体統計
        stats['total_stats'].append({
            'pdb_id': pdb_id,
            'parser': parser,
            'total_bp_count': total_bp,
            'total_canonical_bp_count': total_canonical,
            'total_canonical_ratio': result['total_canonical_ratio'],
            'layer_count': total_layers
        })
        
        # 構造情報
        stats['structure_info'].append({
            'pdb_id': pdb_id,
            'parser': parser,
            'has_pseudoknot': total_layers > 1,
            'pseudoknot_layer_count': max(0, total_layers - 1)  # Layer 0を除く
        })
        
        # 各レイヤーの統計
        layers = result.get('layers', [])
        
        for layer in layers:
            layer_id = layer['layer_id']
            canonical_count = layer['canonical_bp_count']
            total_count = layer['total_bp_count']
            
            layer_info = {
                'pdb_id': pdb_id,
                'parser': parser,
                'layer_id': layer_id,
                'canonical_bp_count': canonical_count,
                'total_bp_count': total_count,
                'canonical_ratio': layer['canonical_ratio']
            }
            
            if layer_id == 0:
                # Main Layer (Core Layer)
                stats['main_layer_stats'].append(layer_info)
            else:
                # Pseudoknot Layers
                stats['pseudoknot_layer_stats'].append(layer_info)
    
    return stats

def plot_main_layer_analysis(stats, output_dir):
    """
    Main Layer (Layer 0) の解析をプロット
    
    Args:
        stats (dict): 統計情報
        output_dir (Path): 出力ディレクトリ
    """
    main_df = pd.DataFrame(stats['main_layer_stats'])
    
    if main_df.empty:
        print("No Main Layer data found!")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Main Layer (Core Layer) Analysis', fontsize=16, fontweight='bold')
    
    # 1. Main Layerのcanonical base pair数の分布
    axes[0, 0].hist(main_df['canonical_bp_count'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 0].set_xlabel('Canonical Base Pairs in Main Layer')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Canonical BP Count in Main Layer')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 統計情報を追加
    mean_canonical = main_df['canonical_bp_count'].mean()
    median_canonical = main_df['canonical_bp_count'].median()
    axes[0, 0].axvline(mean_canonical, color='red', linestyle='--', label=f'Mean: {mean_canonical:.1f}')
    axes[0, 0].axvline(median_canonical, color='green', linestyle='--', label=f'Median: {median_canonical:.1f}')
    axes[0, 0].legend()
    
    # 2. パーサー別比較（Main Layer）
    if 'parser' in main_df.columns and len(main_df['parser'].unique()) > 1:
        parser_comparison = main_df.groupby('parser')['canonical_bp_count'].agg(['mean', 'std', 'count'])
        parser_comparison.plot(kind='bar', y='mean', yerr='std', ax=axes[0, 1], 
                              color=['lightcoral', 'lightblue'], alpha=0.8)
        axes[0, 1].set_title('Parser Comparison: Main Layer Canonical BP Count')
        axes[0, 1].set_ylabel('Mean Canonical BP Count')
        axes[0, 1].tick_params(axis='x', rotation=45)
        axes[0, 1].grid(True, alpha=0.3)
    else:
        axes[0, 1].text(0.5, 0.5, 'Single Parser Data', ha='center', va='center', transform=axes[0, 1].transAxes)
        axes[0, 1].set_title('Parser Comparison: Not Available')
    
    # 3. Main Layer BP数 vs 全体BP数の散布図
    total_df = pd.DataFrame(stats['total_stats'])
    merged_df = main_df.merge(total_df, on=['pdb_id', 'parser'], how='inner')
    
    axes[1, 0].scatter(merged_df['total_bp_count'], merged_df['canonical_bp_count'], 
                      alpha=0.6, color='purple')
    axes[1, 0].set_xlabel('Total Base Pairs in Structure')
    axes[1, 0].set_ylabel('Canonical Base Pairs in Main Layer')
    axes[1, 0].set_title('Main Layer vs Total Structure Size')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 相関係数を計算
    correlation = merged_df['total_bp_count'].corr(merged_df['canonical_bp_count'])
    axes[1, 0].text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                   transform=axes[1, 0].transAxes, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # 4. Main Layer比率（Main Layer BP / Total BP）
    merged_df['main_layer_ratio'] = merged_df['canonical_bp_count'] / merged_df['total_bp_count']
    
    axes[1, 1].hist(merged_df['main_layer_ratio'], bins=25, alpha=0.7, color='orange', edgecolor='black')
    axes[1, 1].set_xlabel('Main Layer Ratio (Main Layer BP / Total BP)')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Distribution of Main Layer Ratio')
    axes[1, 1].grid(True, alpha=0.3)
    
    mean_ratio = merged_df['main_layer_ratio'].mean()
    axes[1, 1].axvline(mean_ratio, color='red', linestyle='--', label=f'Mean: {mean_ratio:.3f}')
    axes[1, 1].legend()
    
    plt.tight_layout()
    output_file = output_dir / 'main_layer_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Main layer analysis plot saved to: {output_file}")

def plot_pseudoknot_layer_analysis(stats, output_dir):
    """
    Pseudoknot Layers (Layer 1+) の解析をプロット
    
    Args:
        stats (dict): 統計情報
        output_dir (Path): 出力ディレクトリ
    """
    pk_df = pd.DataFrame(stats['pseudoknot_layer_stats'])
    struct_df = pd.DataFrame(stats['structure_info'])
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Pseudoknot Layer Analysis', fontsize=16, fontweight='bold')
    
    # 1. Pseudoknotを持つ構造の割合
    has_pk_count = struct_df['has_pseudoknot'].sum()
    total_count = len(struct_df)
    no_pk_count = total_count - has_pk_count
    
    pk_ratio_data = [no_pk_count, has_pk_count]
    pk_ratio_labels = ['No Pseudoknot', 'Has Pseudoknot']
    colors = ['lightgray', 'salmon']
    
    axes[0, 0].pie(pk_ratio_data, labels=pk_ratio_labels, colors=colors, autopct='%1.1f%%', startangle=90)
    axes[0, 0].set_title(f'Pseudoknot Presence\n(Total: {total_count} structures)')
    
    # 2. Pseudoknot layer数の分布
    pk_counts = struct_df[struct_df['has_pseudoknot']]['pseudoknot_layer_count']
    if not pk_counts.empty:
        axes[0, 1].hist(pk_counts, bins=range(int(pk_counts.max()) + 2), alpha=0.7, 
                       color='lightgreen', edgecolor='black', align='left')
        axes[0, 1].set_xlabel('Number of Pseudoknot Layers')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Distribution of Pseudoknot Layer Count')
        axes[0, 1].grid(True, alpha=0.3)
    else:
        axes[0, 1].text(0.5, 0.5, 'No Pseudoknot Structures', ha='center', va='center', transform=axes[0, 1].transAxes)
    
    # 3. 各Pseudoknot layerのcanonical BP数
    if not pk_df.empty:
        layer_bp_counts = pk_df.groupby('layer_id')['canonical_bp_count'].agg(['mean', 'std', 'count'])
        
        x_pos = range(len(layer_bp_counts))
        means = layer_bp_counts['mean']
        stds = layer_bp_counts['std'].fillna(0)
        
        axes[1, 0].bar(x_pos, means, yerr=stds, alpha=0.7, color='lightblue', 
                      capsize=5, error_kw={'linewidth': 2})
        axes[1, 0].set_xlabel('Pseudoknot Layer ID')
        axes[1, 0].set_ylabel('Mean Canonical BP Count')
        axes[1, 0].set_title('Average Canonical BP Count per Pseudoknot Layer')
        axes[1, 0].set_xticks(x_pos)
        axes[1, 0].set_xticklabels([f'Layer {i}' for i in layer_bp_counts.index])
        axes[1, 0].grid(True, alpha=0.3)
        
        # 各レイヤーの構造数を表示
        for i, (idx, row) in enumerate(layer_bp_counts.iterrows()):
            axes[1, 0].text(i, means.iloc[i] + stds.iloc[i] + 0.5, f'n={int(row["count"])}', 
                           ha='center', va='bottom', fontsize=10)
    else:
        axes[1, 0].text(0.5, 0.5, 'No Pseudoknot Layer Data', ha='center', va='center', transform=axes[1, 0].transAxes)
    
    # 4. Pseudoknot層の深さ vs canonical BP数の関係
    if not pk_df.empty:
        scatter_data = pk_df.groupby(['pdb_id', 'layer_id']).first().reset_index()
        
        axes[1, 1].scatter(scatter_data['layer_id'], scatter_data['canonical_bp_count'], 
                          alpha=0.6, color='purple')
        axes[1, 1].set_xlabel('Pseudoknot Layer ID')
        axes[1, 1].set_ylabel('Canonical BP Count')
        axes[1, 1].set_title('Pseudoknot Layer Depth vs BP Count')
        axes[1, 1].grid(True, alpha=0.3)
        
        # トレンドラインを追加
        if len(scatter_data) > 1:
            z = np.polyfit(scatter_data['layer_id'], scatter_data['canonical_bp_count'], 1)
            p = np.poly1d(z)
            axes[1, 1].plot(scatter_data['layer_id'], p(scatter_data['layer_id']), "r--", alpha=0.8)
    else:
        axes[1, 1].text(0.5, 0.5, 'No Pseudoknot Layer Data', ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    output_file = output_dir / 'pseudoknot_layer_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Pseudoknot layer analysis plot saved to: {output_file}")

def plot_parser_comparison(rnaview_stats, dssr_stats, output_dir):
    """
    RNAViewとDSSRの比較プロット
    
    Args:
        rnaview_stats (dict): RNAView統計
        dssr_stats (dict): DSSR統計
        output_dir (Path): 出力ディレクトリ
    """
    # データを結合
    all_stats = {
        'main_layer_stats': rnaview_stats['main_layer_stats'] + dssr_stats['main_layer_stats'],
        'total_stats': rnaview_stats['total_stats'] + dssr_stats['total_stats'],
        'structure_info': rnaview_stats['structure_info'] + dssr_stats['structure_info']
    }
    
    main_df = pd.DataFrame(all_stats['main_layer_stats'])
    total_df = pd.DataFrame(all_stats['total_stats'])
    struct_df = pd.DataFrame(all_stats['structure_info'])
    
    if main_df.empty:
        print("No data available for parser comparison!")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('RNAView vs DSSR Comparison', fontsize=16, fontweight='bold')
    
    # 1. パーサー別Main Layer canonical BP数の比較
    main_df.boxplot(column='canonical_bp_count', by='parser', ax=axes[0, 0])
    axes[0, 0].set_title('Main Layer Canonical BP Count')
    axes[0, 0].set_xlabel('Parser')
    axes[0, 0].set_ylabel('Canonical BP Count')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. パーサー別全体canonical ratio比較
    total_df.boxplot(column='total_canonical_ratio', by='parser', ax=axes[0, 1])
    axes[0, 1].set_title('Total Canonical Ratio')
    axes[0, 1].set_xlabel('Parser')
    axes[0, 1].set_ylabel('Canonical Ratio')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. パーサー別pseudoknot検出率
    pk_detection = struct_df.groupby('parser')['has_pseudoknot'].agg(['sum', 'count'])
    pk_detection['ratio'] = pk_detection['sum'] / pk_detection['count']
    
    pk_detection['ratio'].plot(kind='bar', ax=axes[1, 0], color=['lightcoral', 'lightblue'], alpha=0.8)
    axes[1, 0].set_title('Pseudoknot Detection Rate')
    axes[1, 0].set_ylabel('Detection Rate')
    axes[1, 0].tick_params(axis='x', rotation=45)
    axes[1, 0].grid(True, alpha=0.3)
    
    # 検出率の数値を表示
    for i, (parser, ratio) in enumerate(pk_detection['ratio'].items()):
        axes[1, 0].text(i, ratio + 0.01, f'{ratio:.3f}', ha='center', va='bottom')
    
    # 4. パーサー別layer数の比較
    total_df.boxplot(column='layer_count', by='parser', ax=axes[1, 1])
    axes[1, 1].set_title('Pseudoknot Layer Count')
    axes[1, 1].set_xlabel('Parser')
    axes[1, 1].set_ylabel('Layer Count')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = output_dir / 'parser_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Parser comparison plot saved to: {output_file}")

def generate_summary_statistics(stats, parser_name):
    """
    統計サマリーを生成
    
    Args:
        stats (dict): 統計情報
        parser_name (str): パーサー名
    """
    print(f"\n=== {parser_name.upper()} Summary Statistics ===")
    
    # 全体統計
    total_df = pd.DataFrame(stats['total_stats'])
    if not total_df.empty:
        print(f"Total structures analyzed: {len(total_df)}")
        print(f"Average total BP per structure: {total_df['total_bp_count'].mean():.1f}")
        print(f"Average canonical ratio: {total_df['total_canonical_ratio'].mean():.3f}")
        print(f"Average layer count: {total_df['layer_count'].mean():.1f}")
    
    # Main Layer統計
    main_df = pd.DataFrame(stats['main_layer_stats'])
    if not main_df.empty:
        print(f"\nMain Layer Statistics:")
        print(f"Average canonical BP in main layer: {main_df['canonical_bp_count'].mean():.1f}")
        print(f"Main layer canonical BP range: {main_df['canonical_bp_count'].min():.0f} - {main_df['canonical_bp_count'].max():.0f}")
    
    # Pseudoknot統計
    struct_df = pd.DataFrame(stats['structure_info'])
    if not struct_df.empty:
        pk_structures = struct_df['has_pseudoknot'].sum()
        pk_ratio = pk_structures / len(struct_df)
        print(f"\nPseudoknot Statistics:")
        print(f"Structures with pseudoknots: {pk_structures}/{len(struct_df)} ({pk_ratio:.1%})")
        
        if pk_structures > 0:
            pk_layers = struct_df[struct_df['has_pseudoknot']]['pseudoknot_layer_count']
            print(f"Average pseudoknot layers (when present): {pk_layers.mean():.1f}")
            print(f"Max pseudoknot layers: {pk_layers.max():.0f}")

def main():
    """メイン実行関数"""
    print("=== Pseudoknot Layer Plotting ===")
    
    # 出力ディレクトリを作成
    output_dir = analysis_dir / "plots"
    output_dir.mkdir(exist_ok=True)
    
    # データを読み込み
    rnaview_results = load_analysis_results("rnaview")
    dssr_results = load_analysis_results("dssr")
    
    if not rnaview_results and not dssr_results:
        print("No analysis results found! Please run pseudoknotlayer_analysis.py first.")
        return
    
    # 統計を抽出
    if rnaview_results:
        rnaview_stats = extract_layer_statistics(rnaview_results)
        generate_summary_statistics(rnaview_stats, "rnaview")
        
        # RNAViewプロット
        print("\nGenerating RNAView plots...")
        plot_main_layer_analysis(rnaview_stats, output_dir)
        plot_pseudoknot_layer_analysis(rnaview_stats, output_dir)
    
    if dssr_results:
        dssr_stats = extract_layer_statistics(dssr_results)
        generate_summary_statistics(dssr_stats, "dssr")
        
        # DSSRプロット
        print("\nGenerating DSSR plots...")
        plot_main_layer_analysis(dssr_stats, output_dir)
        plot_pseudoknot_layer_analysis(dssr_stats, output_dir)
    
    # 比較プロット
    if rnaview_results and dssr_results:
        print("\nGenerating parser comparison plots...")
        plot_parser_comparison(rnaview_stats, dssr_stats, output_dir)
    
    print(f"\nAll plots saved to: {output_dir}")
    print("Analysis complete!")

if __name__ == "__main__":
    main()
