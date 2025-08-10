import os
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pathlib import Path
from collections import defaultdict
from scipy.stats import mannwhitneyu
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# 設定
parsers = ['rnaview', 'dssr']
base_dir = Path('analysis')
output_dir = Path('analysis/graphs/0810')
os.makedirs(output_dir, exist_ok=True)  # 出力ディレクトリを作成
os.makedirs(output_dir / 'rnaview', exist_ok=True)  # rnaview 用のサブディレクトリ
os.makedirs(output_dir / 'dssr', exist_ok=True)  # dssr

# フォント設定の改善
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

plt.rcParams.update({
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


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


def plot_non_canonical_ratio_box(
    df: pd.DataFrame,
    parser: str,
    variant: str,
    output_dir: Path = Path("analysis/graphs/0810"),
    seed: int = 0,
    add_stats: bool = True,
):
    """
    Core / Pseudoknot の非カノニカルBP比率の箱ひげ図（モノクロ、論文向け）。
    - 平均は小さな点で表示（凡例なし）
    - y=0..1（統計注記ありのときは少しだけ頭上に余白）
    - （任意）Mann–Whitney U の p値 と Cliff's δ を上部に注記
    """
    # ---- データ整形 ----
    core_total = df['canonical_bp_in_main_layer'] + df['non_canonical_bp_in_main_layer']
    pk_total   = df['canonical_bp_in_pk_layer']  + df['non_canonical_bp_in_pk_layer']

    core_mask = core_total > 0
    core_ratio = (df.loc[core_mask, 'non_canonical_bp_in_main_layer'] / core_total[core_mask]).dropna()

    pk_mask = (df['num_of_layers'] > 1) & (pk_total > 0)
    pk_ratio = (df.loc[pk_mask, 'non_canonical_bp_in_pk_layer'] / pk_total[pk_mask]).dropna()

    if core_ratio.empty and pk_ratio.empty:
        print(f"[plot_non_canonical_ratio_box] 有効データなし ({parser}, {variant})")
        return None

    series_list, xticklabels = [], []
    if not core_ratio.empty:
        series_list.append(core_ratio)
        xticklabels.append(f"Core (n={len(core_ratio)})")
    if not pk_ratio.empty:
        series_list.append(pk_ratio)
        xticklabels.append(f"Pseudoknot (n={len(pk_ratio)})")

    # ---- 作図 ----
    fig, ax = plt.subplots(figsize=(3.35, 3.1))  # 1カラム幅想定
    ax.grid(axis='y', linestyle=(0, (2, 3)), linewidth=0.8, alpha=0.45)
    ax.set_axisbelow(True)

    box = ax.boxplot(
        series_list,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops={'color': '0.1', 'linewidth': 1.6},
    )
    fills = ['0.80', '0.65'][:len(series_list)]
    for b, fc in zip(box['boxes'], fills):
        b.set(facecolor=fc, edgecolor='0.25', linewidth=1.2)

    # 平均値（mean）は視認性を高めるため黒のダイヤマーカーで表示 (凡例付き)
    means = [s.mean() for s in series_list]
    for i, m in enumerate(means, start=1):
        ax.plot(i, m, marker='D', markersize=6, color='black', markeredgecolor='white', markeredgewidth=0.5, label='Mean' if i == 1 else None, zorder=4)
        ax.text(i + 0.18, m, f"{m:.2f}", va='center', ha='left', fontsize=7, color='0.25')

    # 既存の散布図（個々のデータ点）
    rng = np.random.default_rng(seed)
    for i, s in enumerate(series_list, start=1):
        x = np.full(len(s), i, dtype=float) + rng.uniform(-0.06, 0.06, size=len(s))
        ax.scatter(x, s.values, s=6, c='0.2', alpha=0.12, linewidth=0, zorder=1)

    # 凡例（Mean のみ）
    if len(series_list) > 0:
        ax.legend(frameon=False, fontsize=8, loc='upper left')

    ax.set_xticks(range(1, len(series_list) + 1))
    ax.set_xticklabels(xticklabels)
    ax.set_ylabel('Non-canonical base-pair ratio')

    # ▼ タイトルと注記がかぶらないよう headroom を確保
    y_upper = 1.0
    if add_stats and len(series_list) == 2:
        y_upper = 1.06  # ちょい上に余白
    ax.set_ylim(0, y_upper)
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_title('Non-canonical base-pair ratio by layer', fontsize=11, pad=6)

    # ---- 統計（任意） ----
    if add_stats and mannwhitneyu is not None and len(series_list) == 2:
        a, b = series_list
        if len(a) > 0 and len(b) > 0:
            u, p = mannwhitneyu(a, b, alternative='two-sided')
            n, m = len(a), len(b)
            delta = 2 * u / (n * m) - 1  # Cliff's δ

            # 図の最上部から一定マージンを切って固定配置（タイトルと衝突しない）
            y_top = y_upper - 0.035
            ax.plot([1, 1, 2, 2],
                    [y_top - 0.015, y_top, y_top, y_top - 0.015],
                    c='0.2', lw=1)
            ax.text(1.5, y_top + 0.004, f"MWU p={p:.2e}, δ={delta:.2f}",
                    ha='center', va='bottom', fontsize=9)

    # タイトル分の余白を確保
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))

    # ---- 保存：PNG ----
    outdir = (output_dir / parser)
    outdir.mkdir(parents=True, exist_ok=True)
    base = outdir / f"box_noncanonical_ratio_core_vs_pseudoknot_{variant}"
    fig.savefig(str(base.with_suffix(".png")), dpi=600, bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f"Saved PNG: {base.with_suffix('.png')}")
    return base.with_suffix(".png")

# 出力先ディレクトリ作成
for parser in parsers:
    (output_dir / parser).mkdir(parents=True, exist_ok=True)

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

            for multilayer in [True, False]:
                if multilayer:
                    sub = df[df['is_multilayer'] == multilayer]
                else:
                    sub = df
                can_bp_mainlayer = sub['canonical_bp_in_main_layer'] > 0
                non_can_bp_mainlayer = sub['non_canonical_bp_in_main_layer'] > 0
                mainlayer_counts = pd.Series({
                    'Canonical': can_bp_mainlayer.sum(),
                    'Non-Canonical': non_can_bp_mainlayer.sum()
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
                    # pandas の型チェッカ回避のため plt.pie を直接使用
                    values = counts.values
                    plt.pie(
                        values,
                        labels=['Canonical', 'Non-Canonical'],
                        autopct='%1.1f%%',
                        startangle=90
                    )
                    plt.title(f'{title} {parser} {"Multilayer" if multilayer else "Single-layer"} ')
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
        # 追加呼び出し: non-canonical 比率箱ひげ図
        plot_non_canonical_ratio_box(df, parser, variant)