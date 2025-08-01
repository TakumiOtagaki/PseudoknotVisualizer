import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pathlib import Path
from collections import defaultdict

# 設定
parsers = ['rnaview', 'dssr']
base_dir = Path('analysis')
output_dir = Path('analysis/graphs')

# フォント設定の改善
def setup_fonts():
    """利用可能なフォントを確認して設定"""
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    
    # 優先フォント候補
    preferred_fonts = ['Times New Roman', 'Times', 'Liberation Serif', 'DejaVu Serif', 'serif']
    
    for font in preferred_fonts:
        if font in available_fonts:
            plt.rcParams['font.family'] = font
            print(f"Using font: {font}")
            break
    else:
        # フォールバック: デフォルトのserifフォント
        plt.rcParams['font.family'] = 'serif'
        print("Using default serif font")

# フォント設定を実行
setup_fonts()

# matplotlib configs
#フォント設定
#plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
plt.rcParams["font.size"] = 15 # 全体のフォントサイズが変更されます。
#plt.rcParams['xtick.labelsize'] = 9 # 軸だけ変更されます。
#plt.rcParams['ytick.labelsize'] = 24 # 軸だけ変更されます


#軸設定
plt.rcParams['xtick.direction'] = 'in' #x軸の目盛りの向き
plt.rcParams['ytick.direction'] = 'in' #y軸の目盛りの向き
#plt.rcParams['axes.grid'] = True # グリッドの作成
#plt.rcParams['grid.linestyle']='--' #グリッドの線種
plt.rcParams["xtick.minor.visible"] = True  #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True  #y軸補助目盛りの追加
plt.rcParams['xtick.top'] = True  #x軸の上部目盛り
plt.rcParams['ytick.right'] = True  #y軸の右部目盛り


#軸大きさ
#plt.rcParams["xtick.major.width"] = 1.0             #x軸主目盛り線の線幅
#plt.rcParams["ytick.major.width"] = 1.0             #y軸主目盛り線の線幅
#plt.rcParams["xtick.minor.width"] = 1.0             #x軸補助目盛り線の線幅
#plt.rcParams["ytick.minor.width"] = 1.0             #y軸補助目盛り線の線幅
#plt.rcParams["xtick.major.size"] = 10               #x軸主目盛り線の長さ
#plt.rcParams["ytick.major.size"] = 10               #y軸主目盛り線の長さ
#plt.rcParams["xtick.minor.size"] = 5                #x軸補助目盛り線の長さ
#plt.rcParams["ytick.minor.size"] = 5                #y軸補助目盛り線の長さ
#plt.rcParams["axes.linewidth"] = 1.0                #囲みの太さ


#凡例設定
plt.rcParams["legend.fancybox"] = False  # 丸角OFF
plt.rcParams["legend.framealpha"] = 1  # 透明度の指定、0で塗りつぶしなし
plt.rcParams["legend.edgecolor"] = 'black'  # edgeの色を変更
plt.rcParams["legend.markerscale"] = 5 #markerサイズの倍率


# 出力先ディレクトリ作成
for parser in parsers:
    (output_dir / parser).mkdir(parents=True, exist_ok=True)

def create_noncanonical_ratio_boxplot(df, parser, variant, multilayer_only=False):
    """
    各エントリについてCore layerとPseudoknot layerにおけるnon-canonicalの割合を計算してbox plotを作成
    
    Args:
        df: DataFrame with the analysis data
        parser: Parser name (rnaview or dssr)
        variant: Variant name (canonical_only or all)
        multilayer_only: If True, include only multilayer structures; if False, include all structures
    """
    # データフィルタリング
    if multilayer_only:
        sub_df = df[df['is_multilayer'] == True].copy()
        layer_type = "Multilayer"
    else:
        sub_df = df.copy()
        layer_type = "All structures"
    
    if len(sub_df) == 0:
        print(f"No data available for {parser} ({variant}, {layer_type})")
        return
    
    # Main layerでのnon-canonical割合を計算
    main_layer_ratios = []
    pk_layer_ratios = []
    
    for idx, row in sub_df.iterrows():
        # Main layer (Core layer)の割合
        total_main = row['canonical_bp_in_main_layer'] + row['non_canonical_bp_in_main_layer']
        if total_main > 0:
            main_ratio = row['non_canonical_bp_in_main_layer'] / total_main
            main_layer_ratios.append(main_ratio)
        
        # Pseudoknot layerの割合（multilayer構造のみ）
        if multilayer_only or row['num_of_layers'] > 1:
            total_pk = row['canonical_bp_in_pk_layer'] + row['non_canonical_bp_in_pk_layer']
            if total_pk > 0:
                pk_ratio = row['non_canonical_bp_in_pk_layer'] / total_pk
                pk_layer_ratios.append(pk_ratio)
    
    # Box plotのデータ準備
    plot_data = []
    labels = []
    
    if main_layer_ratios:
        plot_data.append(main_layer_ratios)
        labels.append('Core Layer')
    
    if pk_layer_ratios:
        plot_data.append(pk_layer_ratios)
        labels.append('Pseudoknot Layer')
    
    if not plot_data:
        print(f"No valid data for box plot: {parser} ({variant}, {layer_type})")
        return
    
    # Box plot作成
    plt.figure(figsize=(8, 6))
    box_plot = plt.boxplot(plot_data, labels=labels, patch_artist=True)
    
    # 色の設定
    colors = ['lightblue', 'lightcoral']
    for patch, color in zip(box_plot['boxes'], colors[:len(plot_data)]):
        patch.set_facecolor(color)
    
    # 中央線（median）の色を黒に設定
    for median in box_plot['medians']:
        median.set_color('black')
        median.set_linewidth(2)
    
    # 平均値を計算して表示
    for i, data in enumerate(plot_data):
        mean_val = np.mean(data)
        plt.scatter(i + 1, mean_val, color='red', marker='D', s=15, zorder=5, label='Mean' if i == 0 else "")
        plt.text(i + 1, mean_val + 0.02, f'{mean_val:.3f}', ha='center', va='bottom', fontweight='bold', color='red')
    
    # 凡例を追加（平均値用）
    plt.legend(loc='upper right')
    
    plt.ylabel('Non-canonical Base Pair Ratio')
    
    # タイトルにエントリ数を追加
    entry_counts = [len(data) for data in plot_data]
    entry_info = ", ".join([f"{labels[i]}: {count}" for i, count in enumerate(entry_counts)])
    plt.title(f'{parser.upper()} Non-canonical BP Ratio ({variant}, {layer_type})\nEntries - {entry_info}')
    
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # ファイル名生成と保存
    multilayer_suffix = "multilayer" if multilayer_only else "including_singlelayer"
    fn = f'boxplot_noncanonical_ratio_{multilayer_suffix}_{variant}.png'
    plt.savefig(os.path.join(output_dir, parser, fn))
    print(f"Saved box plot for {parser} ({variant}, {layer_type}) to {output_dir / parser / fn}")
    print(f"  Core Layer: {len(main_layer_ratios)} structures")
    print(f"  Pseudoknot Layer: {len(pk_layer_ratios)} structures")
    plt.close()

def create_weighted_pie_charts(df, parser, variant, multilayer_only=False):
    """
    塩基対数ベースの重み付けでpie chartを作成
    （構造の大きさに比例した重み付け）
    
    Args:
        df: DataFrame with the analysis data
        parser: Parser name (rnaview or dssr)
        variant: Variant name (canonical_only or all)
        multilayer_only: If True, include only multilayer structures; if False, include all structures
    """
    # データフィルタリング
    if multilayer_only:
        sub_df = df[df['is_multilayer'] == True].copy()
        layer_type = "Multilayer"
    else:
        sub_df = df.copy()
        layer_type = "All structures"
    
    if len(sub_df) == 0:
        print(f"No data available for weighted pie chart: {parser} ({variant}, {layer_type})")
        return
    
    # 塩基対数ベースで重み付けして計算
    mainlayer_counts = pd.Series({
        'Canonical': sub_df['canonical_bp_in_main_layer'].sum(),
        'Non-Canonical': sub_df['non_canonical_bp_in_main_layer'].sum()
    })
    
    pseudoknotlayer_counts = pd.Series({
        'Canonical': sub_df['canonical_bp_in_pk_layer'].sum(),
        'Non-Canonical': sub_df['non_canonical_bp_in_pk_layer'].sum()
    })
    
    # 各レイヤーのpie chartを作成
    for counts, title in zip(
        [mainlayer_counts, pseudoknotlayer_counts],
        ['MainLayer', 'PseudoknotLayer']
    ):
        # データが存在する場合のみプロット
        if counts.sum() > 0:
            plt.figure()
            counts.plot.pie(
                autopct='%1.1f%%',
                labels=['Canonical', 'Non-Canonical'],
                startangle=90,
                title=f'{title} {parser.upper()} ({layer_type}) - Weighted by BP count\nTotal BPs: {int(counts.sum())}'
            )
            plt.ylabel('')
            plt.tight_layout()
            
            multilayer_suffix = "multilayer" if multilayer_only else "including_singlelayer"
            fn = f'pie_weighted_{title}_{multilayer_suffix}_{variant}.png'
            plt.savefig(os.path.join(output_dir, parser, fn))
            print(f"Saved weighted pie chart for {parser} ({variant}, {layer_type}, {title}) to {output_dir / parser / fn}")
            print(f"  Total BPs: {int(counts.sum())}, Canonical: {counts['Canonical']}, Non-Canonical: {counts['Non-Canonical']}")
            plt.close()
        else:
            print(f"No base pairs found for {title} in {parser} ({variant}, {layer_type})")


for parser in parsers:
    for variant in ['canonical_only', 'all']:

        total_bps, total_canonical_bps, total_noncanonical_bps = 0, 0, 0

        # JSON読み込み
        path = base_dir / f'pseudoknot_analysis_{parser}_{variant}.json'
        with open(path) as f:
            data = json.load(f)
        # layers --> pdb_id v.s. layer_id, canonical_bp_count, non_canonical_bp_count
        data_list = list()
        for item in data: # 各 pdb ごとに分解していく
            # sys.exit()
            pdb_id = item['pdb_id']
            layers = item['layers']
            num_of_layers = item['pseudoknot_layer_count']
            if num_of_layers == 0:
                # print(f"Warning: {pdb_id} has no pseudoknot layers.")
                num_of_layers = 1  # レイヤーがない場合は1とみなす
                continue  # レイヤーがない場合はスキップ...
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
            data_list.append(
                {
                    "pdb_id": pdb_id,
                    "num_of_layers": num_of_layers,
                    "total_bp_count": total_bp_count,
                    "canonical_bp_in_main_layer": canonical_bp_in_main_layer,
                    "non_canonical_bp_in_main_layer": non_canonical_bp_in_main_layer,
                    "canonical_bp_in_pk_layer": canonical_bp_in_pk_layer,
                    "non_canonical_bp_in_pk_layer": non_canonical_bp_in_pk_layer,
                }
            )
            total_bps += total_bp_count
            total_canonical_bps += canonical_bp_in_main_layer + canonical_bp_in_pk_layer
            total_noncanonical_bps += non_canonical_bp_in_main_layer + non_canonical_bp_in_pk_layer
        print(f"Total BPs: {total_bps}, Total Canonical BPs: {total_canonical_bps}, Total Non-Canonical BPs: {total_noncanonical_bps}, for {parser} ({variant})")
        df = pd.DataFrame(data_list)
        df = df.set_index("pdb_id")
        # print(f"DataFrame for {parser} ({variant}):\n", df.head())

        # (a) トップレイヤーにcanonicalがあるか の円グラフ
        if variant == "all":
            df['is_multilayer'] = df["num_of_layers"] > 1 # 元々 multilayer かどうかの場合分けで、2 通り plotting
            print(f"number of multilayer structures: {df['is_multilayer'].sum()} out of {len(df)} total structures, for {parser} ({variant})")

            # Box plot for non-canonical ratios
            create_noncanonical_ratio_boxplot(df, parser, variant, multilayer_only=True)
            create_noncanonical_ratio_boxplot(df, parser, variant, multilayer_only=False)

            # Weighted pie charts (新しく追加)
            create_weighted_pie_charts(df, parser, variant, multilayer_only=True)
            create_weighted_pie_charts(df, parser, variant, multilayer_only=False)

            for multilayer in [True, False]:
                if multilayer:
                    sub = df[df['is_multilayer'] == multilayer]
                else:
                    sub = df
                mainlayer_counts = pd.Series({
                    'Canonical': sub['canonical_bp_in_main_layer'].sum(),
                    'Non-Canonical': sub['non_canonical_bp_in_main_layer'].sum()
                })
                pseudoknotlayer_counts = pd.Series({
                    'Canonical': sub['canonical_bp_in_pk_layer'].sum(),
                    'Non-Canonical': sub['non_canonical_bp_in_pk_layer'].sum()
                })

                for counts, title in zip(
                    [mainlayer_counts, pseudoknotlayer_counts],
                    ['MainLayer', 'PseudoknotLayer']
                ):
                    plt.figure()
                    counts.plot.pie(
                        autopct='%1.1f%%',
                        labels=['Canonical', 'Non-Canonical'],
                        startangle=90,
                        title=f'{title} {parser} {"Multilayer" if multilayer else "Single-layer"} '
                    )
                    plt.ylabel('')
                    plt.tight_layout()
                    fn = f'pie_{title}{"multilayer" if multilayer else "including_singlelayer"}_{variant}.png'
                    plt.savefig(os.path.join(output_dir, parser, fn))
                    print(f"Saved pie chart for {parser} ({variant}, {'multilayer' if multilayer else 'single-layer'}) to {output_dir / parser / fn}")
                    plt.close()
        else: 
            pass

        # (b) pseudoknot layer 数の頻度棒グラフ
        freq = df['num_of_layers'].value_counts(normalize=True).sort_index()
        plt.figure()
        ax = freq.plot.bar( color = "orange")
        ax.set_xlabel('Pseudoknot Layer Count')
        ax.set_xticks([int(x) for x in range(freq.index.max() + 1)])
        ax.set_xticklabels([int(x) for x in (range(freq.index.max() + 1))], rotation=0  )
        
        ax.set_xlim(-1, freq.index.max()+1)  # x 軸の範囲を 0 から最大値に設定
        # plt.xlim(0, 10)
        plt.ylabel('Frequency')
        # y の表示領域を 0 から 0.7 に固定
        plt.ylim(0, 0.7)
        plt.title(f'{parser} Pseudoknot Layers ({variant})')
        plt.tight_layout()
        fn = f'bar_pseudoknot_layers_{variant}.png'
        plt.savefig(os.path.join(output_dir, parser, fn))
        print(f"Saved bar graph for {parser} ({variant}) to {output_dir / parser / fn}")
        plt.close()