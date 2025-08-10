#!/usr/bin/env python3
"""
Multi Base Pairing Entries Plotting Script

"(i, j) と (i, j') といった multi base pairing を持ったエントリ"専用の可視化スクリプト

前回のグラフとの差分を見るために、pseudoknotlayer_plotting.pyと同じスタイルで可視化
"""

import os
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pathlib import Path
from collections import defaultdict

def setup_fonts():
    """利用可能なフォントを確認して設定"""
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    preferred_fonts = ['Times New Roman', 'Times', 'Liberation Serif', 'DejaVu Serif', 'serif']
    
    for font in preferred_fonts:
        if font in available_fonts:
            plt.rcParams['font.family'] = font
            print(f"Using font: {font}")
            break
    else:
        plt.rcParams['font.family'] = 'serif'
        print("Using default serif font")

def setup_matplotlib():
    """matplotlibの設定"""
    setup_fonts()
    
    plt.rcParams["font.size"] = 15
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    plt.rcParams["legend.fancybox"] = False
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams["legend.edgecolor"] = 'black'
    plt.rcParams["legend.markerscale"] = 5

def analyze_multi_pairing_data(data):
    """multi base pairingデータを解析してDataFrameに変換"""
    data_list = []
    total_bps, total_canonical_bps, total_noncanonical_bps = 0, 0, 0
    
    for item in data:
        pdb_id = item['pdb_id']
        layers = item['layers']
        num_of_layers = item['pseudoknot_layer_count']
        
        if num_of_layers == 0:
            num_of_layers = 1
            continue
        
        canonical_bp_in_main_layer = 0
        non_canonical_bp_in_main_layer = 0
        canonical_bp_in_pk_layer = 0
        non_canonical_bp_in_pk_layer = 0
        
        for layer in layers:
            layer_id = layer['layer_id']
            if layer_id == 0:  # トップレイヤー
                canonical_bp_in_main_layer = layer['canonical_bp_count']
                non_canonical_bp_in_main_layer = layer['non_canonical_bp_count']
            else:  # PKレイヤー
                canonical_bp_in_pk_layer += layer['canonical_bp_count']
                non_canonical_bp_in_pk_layer += layer['non_canonical_bp_count']
        
        total_bp_count = item['total_bp_count']
        dup_canonical_count = len(item.get('dup_canonical_pairs', []))
        
        data_list.append({
            "pdb_id": pdb_id,
            "num_of_layers": num_of_layers,
            "total_bp_count": total_bp_count,
            "canonical_bp_in_main_layer": canonical_bp_in_main_layer,
            "non_canonical_bp_in_main_layer": non_canonical_bp_in_main_layer,
            "canonical_bp_in_pk_layer": canonical_bp_in_pk_layer,
            "non_canonical_bp_in_pk_layer": non_canonical_bp_in_pk_layer,
            "dup_canonical_count": dup_canonical_count
        })
        
        total_bps += total_bp_count
        total_canonical_bps += canonical_bp_in_main_layer + canonical_bp_in_pk_layer
        total_noncanonical_bps += non_canonical_bp_in_main_layer + non_canonical_bp_in_pk_layer
    
    print(f"Multi Pairing Entries - Total BPs: {total_bps}, Canonical: {total_canonical_bps}, Non-Canonical: {total_noncanonical_bps}")
    
    df = pd.DataFrame(data_list)
    df = df.set_index("pdb_id")
    return df

def create_pie_charts(df, output_dir, title_prefix="MultiPairing"):
    """円グラフを作成"""
    df['is_multilayer'] = df["num_of_layers"] > 1
    print(f"Multi-layer structures: {df['is_multilayer'].sum()} out of {len(df)} total")
    
    for multilayer in [True, False]:
        if multilayer:
            sub = df[df['is_multilayer'] == multilayer]
            suffix = "multilayer"
        else:
            sub = df
            suffix = "including_singlelayer"
        
        if len(sub) == 0:
            print(f"No data for {'multilayer' if multilayer else 'single-layer'} structures")
            continue
        
        # Main Layer円グラフ
        can_bp_mainlayer = sub['canonical_bp_in_main_layer'] > 0
        non_can_bp_mainlayer = sub['non_canonical_bp_in_main_layer'] > 0
        mainlayer_counts = pd.Series({
            'Canonical': can_bp_mainlayer.sum(),
            'Non-Canonical': non_can_bp_mainlayer.sum()
        })
        
        # Pseudoknot Layer円グラフ
        pseudoknotlayer_counts = pd.Series({
            'Canonical': sub['canonical_bp_in_pk_layer'].sum(),
            'Non-Canonical': sub['non_canonical_bp_in_pk_layer'].sum()
        })
        
        for counts, layer_type in zip(
            [mainlayer_counts, pseudoknotlayer_counts],
            ['MainLayer', 'PseudoknotLayer']
        ):
            plt.figure()
            counts.plot.pie(
                autopct='%1.1f%%',
                labels=['Canonical', 'Non-Canonical'],
                startangle=90,
                title=f'{layer_type} {title_prefix} {"Multilayer" if multilayer else "Single-layer"}'
            )
            plt.ylabel('')
            plt.tight_layout()
            
            fn = f'pie_{layer_type}_{suffix}_multipairing.png'
            plt.savefig(output_dir / fn)
            print(f"Saved: {output_dir / fn}")
            plt.close()

def create_layer_distribution(df, output_dir, title_prefix="MultiPairing"):
    """レイヤー数分布の棒グラフを作成"""
    freq = df['num_of_layers'].value_counts(normalize=True).sort_index()
    
    plt.figure()
    ax = freq.plot.bar(color="red", alpha=0.7)  # multi pairingは赤で区別
    ax.set_xlabel('Pseudoknot Layer Count')
    ax.set_xticks([int(x) for x in range(freq.index.max() + 1)])
    ax.set_xticklabels([int(x) for x in range(freq.index.max() + 1)], rotation=0)
    ax.set_xlim(-1, freq.index.max() + 1)
    
    plt.ylabel('Frequency')
    plt.ylim(0, 0.7)
    plt.title(f'{title_prefix} Pseudoknot Layers Distribution')
    plt.tight_layout()
    
    fn = 'bar_pseudoknot_layers_multipairing.png'
    plt.savefig(output_dir / fn)
    print(f"Saved: {output_dir / fn}")
    plt.close()

def create_duplicate_analysis(df, output_dir):
    """重複canonical塩基対の分析グラフ"""
    
    # 重複数の分布
    dup_counts = df['dup_canonical_count'].value_counts().sort_index()
    
    plt.figure()
    ax = dup_counts.plot.bar(color="darkred")
    ax.set_xlabel('Number of Duplicate Canonical Pairs')
    ax.set_ylabel('Number of Entries')
    plt.title('Distribution of Duplicate Canonical Pairs')
    plt.tight_layout()
    
    fn = 'bar_duplicate_canonical_distribution.png'
    plt.savefig(output_dir / fn)
    print(f"Saved: {output_dir / fn}")
    plt.close()
    
    # 重複数 vs レイヤー数の散布図
    plt.figure()
    plt.scatter(df['dup_canonical_count'], df['num_of_layers'], 
                alpha=0.6, color='darkred')
    plt.xlabel('Number of Duplicate Canonical Pairs')
    plt.ylabel('Number of Pseudoknot Layers')
    plt.title('Duplicate Canonical Pairs vs Pseudoknot Layers')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    fn = 'scatter_duplicate_vs_layers.png'
    plt.savefig(output_dir / fn)
    print(f"Saved: {output_dir / fn}")
    plt.close()

def main():
    """メイン実行関数"""
    print("🎨 Multi Base Pairing Entries 可視化スクリプト")
    print("=" * 60)
    
    # 設定
    setup_matplotlib()
    
    # 入力・出力設定
    input_file = Path("analysis/multi_pairing_entries_dssr.json")
    output_dir = Path("analysis/graphs/multi_pairing")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # データ読み込み
    print(f"📖 データ読み込み: {input_file}")
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    print(f"📊 Multi base pairing エントリ数: {len(data)}")
    
    # データ解析
    df = analyze_multi_pairing_data(data)
    
    # 可視化
    print("\n🎨 グラフ作成中...")
    
    # 1. 円グラフ（Main Layer / Pseudoknot Layer別）
    create_pie_charts(df, output_dir, "MultiPairing_DSSR")
    
    # 2. レイヤー数分布
    create_layer_distribution(df, output_dir, "MultiPairing_DSSR")
    
    # 3. 重複canonical塩基対分析
    create_duplicate_analysis(df, output_dir)

    # 4. 箱ひげ図（レイヤー数ごとの塩基対数）
    box_plot_dimensions(df, output_dir)
    
    # 統計サマリー出力
    print("\n📈 統計サマリー")
    print("-" * 40)
    print(f"総エントリ数: {len(df)}")
    print(f"Multi-layerエントリ数: {df['is_multilayer'].sum()}")
    print(f"平均レイヤー数: {df['num_of_layers'].mean():.2f}")
    print(f"平均重複canonical数: {df['dup_canonical_count'].mean():.1f}")
    print(f"平均塩基対数: {df['total_bp_count'].mean():.1f}")
    
    canonical_ratio = (df['canonical_bp_in_main_layer'] + df['canonical_bp_in_pk_layer']).sum() / df['total_bp_count'].sum()
    print(f"Canonical比率: {canonical_ratio:.3f}")
    
    print(f"\n✅ 全てのグラフが {output_dir} に保存されました")

def box_plot_non_canonical_ratio(df, output_dir):
    """レイヤーごと(Main layer or Pseudoknot layer) の非canonical塩基対比率の箱ひげ図を作成"""
    plt.figure(figsize=(10, 6))
    df['non_canonical_ratio_main'] = df['non_canonical_bp_in_main_layer'] / df['total_bp_count']
    df['non_canonical_ratio_pk'] = df['non_canonical_bp_in_pk_layer'] / df['total_bp_count']
    
if __name__ == "__main__":
    main()
