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
from matplotlib.markers import MarkerStyle

# 設定
parsers = ['rnaview', 'dssr']
base_dir = Path('analysis')
output_dir = Path('analysis/graphs/0829')
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
plt.rcParams["font.size"] = 18 # 全体のフォントサイズが変更されます。
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
    output_dir: Path,
    seed: int = 0,
    add_stats: bool = False,
    show_mean: bool = True,
    include_no_pk: bool = False,
):
    """
    Core / Pseudoknot の non-canonical BP 比率の箱ひげ図（論文向けミニマル版）。
    - 右/上スパイン削除、横グリッドのみ
    - 平均は小さな点（●）で表示。凡例や数値注記は出さない
    - （任意）MWU の p値 + Cliff's δ を上部に固定表示
    """
    # ---- データ整形 ----
    core_total = df['canonical_bp_in_main_layer'] + df['non_canonical_bp_in_main_layer']
    pk_total   = df['canonical_bp_in_pk_layer']  + df['non_canonical_bp_in_pk_layer']

    print(f"  [RATIO BOX] {parser} {variant} filtering details:")
    print(f"    - Total structures in input: {len(df)}")
    
    core_mask = core_total > 0
    core_ratio = (df.loc[core_mask, 'non_canonical_bp_in_main_layer'] / core_total[core_mask]).dropna()
    
    print(f"    - Core layer: {core_mask.sum()} structures have >0 BPs, {len(core_ratio)} valid ratios")
    if len(core_ratio) > 0:
        core_canonical = df.loc[core_mask, 'canonical_bp_in_main_layer'].sum()
        core_non_canonical = df.loc[core_mask, 'non_canonical_bp_in_main_layer'].sum()
        print(f"      * Core canonical BPs: {core_canonical}")
        print(f"      * Core non-canonical BPs: {core_non_canonical}")

    pk_mask = (df['num_of_layers'] > 1) & (pk_total > 0)
    pk_ratio = (df.loc[pk_mask, 'non_canonical_bp_in_pk_layer'] / pk_total[pk_mask]).dropna()
    if include_no_pk:
        # PK層がないものも 0 として含める
        no_pk_idx = df.index[pk_total == 0]
        if len(no_pk_idx) > 0:
            zeros = pd.Series(0.0, index=no_pk_idx)
            pk_ratio = pd.concat([pk_ratio, zeros])

    print(f"    - PK layer: {pk_mask.sum()} structures have >1 layers & >0 PK BPs, {len(pk_ratio)} valid ratios")
    if pk_mask.sum() > 0:
        pk_canonical = df.loc[pk_mask, 'canonical_bp_in_pk_layer'].sum()
        pk_non_canonical = df.loc[pk_mask, 'non_canonical_bp_in_pk_layer'].sum()
        print(f"      * PK canonical BPs: {pk_canonical}")
        print(f"      * PK non-canonical BPs: {pk_non_canonical}")

    # 最終的にボックスプロットに使われるチェーン数を表示
    n_core = len(core_ratio)
    n_pk = len(pk_ratio)
    if include_no_pk:
        print(f"    [FINAL] {parser} {variant} (incl. no PK): Box plot will show Core={n_core}, Pseudoknot(incl. none)={n_pk}")
    else:
        print(f"    [FINAL] {parser} {variant}: Box plot will show Core={n_core}, Pseudoknot(≥1 layer)={n_pk}")

    if core_ratio.empty and pk_ratio.empty:
        print(f"    [ERROR] 有効データなし ({parser}, {variant})")
        return None

    series_list, xticklabels = [], []
    if not core_ratio.empty:
        series_list.append(core_ratio)
        # xticklabels.append(f"Core (n={len(core_ratio)})")
        xticklabels.append(f"Core layer")
    if not pk_ratio.empty:
        series_list.append(pk_ratio)
        # xticklabels.append(f"Pseudoknot (n={len(pk_ratio)})")
        xticklabels.append(f"Pseudoknot layer")

    # ---- 作図 ----
    fig, ax = plt.subplots(figsize=(3.35, 3.1), constrained_layout=True)  # 1カラム幅想定
    ax.grid(axis='y', linestyle=(0, (2, 3)), linewidth=0.8, alpha=0.45)
    ax.set_axisbelow(True)

    box = ax.boxplot(
        series_list,
        widths=0.5,
        patch_artist=True,
        showfliers=False,  # 外れ値はジッター点に任せる
        medianprops={'color': '0.1', 'linewidth': 1.0},
        whiskerprops={'linewidth': 1.2, 'color': '0.25'},
        capprops={'linewidth': 1.2, 'color': '0.25'},
    )
    fills = ['0.80', '0.65'][:len(series_list)]
    for b, fc in zip(box['boxes'], fills):
        b.set(facecolor=fc, edgecolor='0.25', linewidth=1.2)

    # 平均（小さな黒点）— 凡例/数値もつける
    # if show_mean:
    #     means = [s.mean() for s in series_list]
    #     ax.scatter(range(1, len(series_list) + 1), means, s=18, marker='*', c='0.1', zorder=3, label='Mean')

    # ax.legend(markerscale=1.5)
    # 平均（小さなマーカー）＋数値ラベル
    if show_mean:
        means = [s.mean() for s in series_list]

        # 平均マーカー（★）
        sc = ax.scatter(range(1, len(series_list) + 1), means,
                        s=24, marker=MarkerStyle('*'), c='0.1', zorder=3, label='Mean')

        # 数値ラベルを右横に
        y0, y1 = ax.get_ylim()
        pad_y = 0.02 * (y1 - y0)  # 上端からの安全マージン
        epsilon = 0.007
        for i, m in enumerate(means, start=1):
            y = min(m + epsilon, y1 - pad_y)  # 文字が上端をはみ出さないように
            ax.annotate(f"{m:.3f}", (i, y),
                        xytext=(6, 0), textcoords="offset points",  # 右に6ptずらす
                        ha='left', va='center', fontsize=8, color='0.25')

        # 凡例（マーカーだけ小さめに）
        ax.legend(handles=[sc],
                  markerscale=0.8, scatterpoints=1,
                  handlelength=0.8, handletextpad=0.4,
                  frameon=False, fontsize=8)

    # ジッター散布（超薄）
    rng = np.random.default_rng(seed)
    for i, s in enumerate(series_list, start=1):
        if len(s) == 0:
            continue
        x = np.full(len(s), i, dtype=float) + rng.uniform(-0.055, 0.055, size=len(s))
        ax.scatter(x, s.values, s=5, c='0.2', alpha=0.10, linewidth=0, zorder=1)

    # 軸・ラベル
    ax.set_xticks(range(1, len(series_list) + 1))
    ax.set_xticklabels(xticklabels)
    ax.set_ylabel('Non-canonical base-pair ratio')
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_title('Non-canonical base-pair ratio by layer', fontsize=14, pad=6)

    # ---- y上端と統計注記（必要なときのみ） ----
    y_upper = 1.0
    if add_stats and mannwhitneyu is not None and len(series_list) == 2:
        a, b = series_list
        if len(a) > 0 and len(b) > 0:
            # Uとp、Cliff's δ
            u, p = mannwhitneyu(a, b, alternative='two-sided')
            n, m = len(a), len(b)
            delta = 2 * u / (n * m) - 1  # U から算出
            y_upper = 1.06  # 注記のために頭上を少し確保
            ax.set_ylim(0, y_upper)

            y_top = y_upper - 0.035  # 固定相対位置（タイトルと衝突しない）
            ax.plot([1, 1, 2, 2], [y_top - 0.015, y_top, y_top, y_top - 0.015],
                    c='0.2', lw=1)
            ax.text(1.5, y_top + 0.004, f"MWU p={p:.2e}, δ={delta:.3f}",
                    ha='center', va='bottom', fontsize=12)
    else:
        ax.set_ylim(0, y_upper)

    # ---- 保存（PNG / 白背景） ----
    outdir = output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    base = outdir / f"box_noncanonical_ratio_core_vs_pseudoknot_{variant}{'_incl_nopk' if include_no_pk else ''}"
    fig.savefig(str(base.with_suffix(".png")), dpi=600, bbox_inches='tight', transparent=False)
    plt.close(fig)
    print(f"Saved PNG: {base.with_suffix('.png')}")
    return base.with_suffix(".png")




def main():

    # 出力先ディレクトリ作成
    for parser in parsers:
        (output_dir / parser).mkdir(parents=True, exist_ok=True)

    for parser in parsers:
        print(f"---- Processing parser: {parser} ----")
        for variant in ['canonical_only', 'all']:

            total_bps, total_canonical_bps, total_noncanonical_bps = 0, 0, 0
            structures_with_layers = 0
            structures_skipped = 0

            # JSON読み込み
            path = base_dir / f'pseudoknot_analysis_{parser}_{variant}.json'
            with open(path) as f:
                data = json.load(f)
            
            print(f"  [FILTER] 初期データ読み込み: {len(data)} structures loaded from JSON file")
            
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
                    structures_skipped += 1
                    continue  # レイヤーがない場合はスキップ...
                
                structures_with_layers += 1
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
            
            print(f"  [FILTER] レイヤーなし構造をスキップ: {structures_skipped} structures skipped (no layers)")
            print(f"  [FILTER] 残存構造数: {structures_with_layers} structures with layers")
            print(f"  [FILTER] 総塩基対数: {total_bps} BPs")
            print(f"  [FILTER] Canonical BPs: {total_canonical_bps} ({total_canonical_bps/total_bps*100:.1f}%)")
            print(f"  [FILTER] Non-canonical BPs: {total_noncanonical_bps} ({total_noncanonical_bps/total_bps*100:.1f}%)")
            df = pd.DataFrame(data_list)
            df = df.set_index("pdb_id")
            
            print(f"  [FILTER] DataFrame作成後: {len(df)} structures in final dataset")
            
            # DataFrame作成後の統計情報
            df_total_canonical = (df['canonical_bp_in_main_layer'] + df['canonical_bp_in_pk_layer']).sum()
            df_total_non_canonical = (df['non_canonical_bp_in_main_layer'] + df['non_canonical_bp_in_pk_layer']).sum()
            df_total_bps = df_total_canonical + df_total_non_canonical
            
            print(f"  [FILTER] DataFrame内統計:")
            print(f"    - Main layer canonical BPs: {df['canonical_bp_in_main_layer'].sum()}")
            print(f"    - Main layer non-canonical BPs: {df['non_canonical_bp_in_main_layer'].sum()}")
            print(f"    - PK layer canonical BPs: {df['canonical_bp_in_pk_layer'].sum()}")
            print(f"    - PK layer non-canonical BPs: {df['non_canonical_bp_in_pk_layer'].sum()}")
            print(f"    - Total canonical BPs: {df_total_canonical} ({df_total_canonical/df_total_bps*100:.1f}%)")
            print(f"    - Total non-canonical BPs: {df_total_non_canonical} ({df_total_non_canonical/df_total_bps*100:.1f}%)")
            # print(f"DataFrame for {parser} ({variant}):\n", df.head())

            # (a) トップレイヤーにcanonicalがあるか の円グラフ
            if variant == "all":
                df['is_multilayer'] = df["num_of_layers"] > 1 # 元々 multilayer かどうかの場合分けで、2 通り plotting
                multilayer_count = df['is_multilayer'].sum()
                total_structures = len(df)
                print(f"  [FILTER] Multilayer分析:")
                print(f"    - Multilayer structures: {multilayer_count} ({multilayer_count/total_structures*100:.1f}%)")
                print(f"    - Single-layer structures: {total_structures - multilayer_count} ({(total_structures - multilayer_count)/total_structures*100:.1f}%)")

                for multilayer in [True, False]:
                    if multilayer:
                        sub = df[df['is_multilayer'] == multilayer]
                        layer_type = "multilayer"
                    else:
                        sub = df
                        layer_type = "all structures"
                    
                    print(f"    [{layer_type.upper()}] Processing {len(sub)} structures:")
                    
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
                    
                    print(f"      Main layer structures with canonical BPs: {mainlayer_counts['Canonical']}")
                    print(f"      Main layer structures with non-canonical BPs: {mainlayer_counts['Non-Canonical']}")
                    print(f"      Total canonical BPs in PK layers: {pseudoknotlayer_counts['Canonical']}")
                    print(f"      Total non-canonical BPs in PK layers: {pseudoknotlayer_counts['Non-Canonical']}")

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
                        parser_text = "DSSR" if parser == "dssr" else "RNAView"
                        plt.title(f'{title} {parser_text} {"Multilayer" if multilayer else "Single-layer"} ')
                        plt.ylabel('')
                        plt.tight_layout()
                        fn = f'pie_{title}{"multilayer" if multilayer else "including_singlelayer"}_{variant}.png'
                        plt.savefig(os.path.join(output_dir, parser, fn))
                        print(f"Saved pie chart for {parser} ({variant}, {'multilayer' if multilayer else 'single-layer'}) to {output_dir / parser / fn}")
                        plt.close()
            else: 
                pass

            # (b) pseudoknot layer 数の頻度棒グラフ + 累積度数ライン
            freq_counts = df['num_of_layers'].value_counts().sort_index()
            freq = (freq_counts / freq_counts.sum()).sort_index()
            cumulative = freq.cumsum()
            
            print(f"  [FILTER] Layer分布統計:")
            for layer_count, count in freq_counts.items():
                percentage = count / freq_counts.sum() * 100
                print(f"    - {layer_count} layers: {count} structures ({percentage:.1f}%)")
            
            plt.figure()
            ax = freq.plot.bar(color="orange", alpha=0.85, width=0.8)
            # 上・右の目盛り線/スパイン不要化
            ax.tick_params(top=False, right=False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(which='both', top=False, right=False)
            ax.minorticks_off()  # minor を完全にオフ
            ax.set_xlabel('Pseudoknot Order (number of layers)', fontsize=16)
            ax.set_xticks([int(x) for x in range(freq.index.max())])
            ax.set_xticklabels([int(x) for x in (range(1, freq.index.max()+1))], rotation=0 , fontsize=15)
            ax.tick_params(axis='y', labelsize=15)
            ax.set_xlim(-0.6, freq.index.max()+0.3)
            plt.ylabel('Frequency', fontsize=16)
            # 累積度数ライン（同じ y 軸: 正規化なので 0→1）
            ax.plot(cumulative.index - 1, cumulative.values, color='gray', marker='o', linewidth=1.2, markersize=5, label='Cumulative', linestyle='--')
            # 各累積点に軽く値を表示（上から重ねて）
            for x, y in zip(cumulative.index, cumulative.values):
                ax.text(x - 1, y + 0.02, f"{y:.3f}", ha='center', va='bottom', fontsize=11, color='black')
            # y の表示領域（最大 1 を超えないように）
            cum_max = float(cumulative.max()) if len(cumulative) else 0.0
            top = max(0.7, cum_max)
            ax.set_ylim(0, min(1.05, top + 0.08))
            parser_text = "DSSR" if parser == "dssr" else "RNAView"
            plt.title(f'{parser_text} Pseudoknot Order ({variant.replace("_", " ")})', fontsize=18)
            ax.legend(frameon=False, fontsize=12, loc='lower right', markerscale=1.0)
            plt.tight_layout()
            fn = f'bar_pseudoknot_layers_{variant}.png'
            plt.savefig(os.path.join(output_dir, parser, fn))
            print(f"Saved bar+cumulative graph for {parser} ({variant}) to {output_dir / parser / fn}")
            plt.close()
            # 追加呼び出し: non-canonical 比率箱ひげ
            # (1) PK>=1 のみ
            plot_non_canonical_ratio_box(df, parser, variant, output_dir=output_dir / parser, show_mean=True, include_no_pk=False)
            # (2) PKなしも0として含める版
            plot_non_canonical_ratio_box(df, parser, variant, output_dir=output_dir / parser, show_mean=True, include_no_pk=True)


if __name__ == "__main__":
    main()