import os
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

# 設定
parsers = ['rnaview', 'dssr']
base_dir = Path('analysis')
output_dir = Path('analysis/graphs')

# 出力先ディレクトリ作成
for parser in parsers:
    (output_dir / parser).mkdir(parents=True, exist_ok=True)

for parser in parsers:
    for variant in ['canonical_only', 'all']:
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
        df = pd.DataFrame(data_list)
        df = df.set_index("pdb_id")
        # print(f"DataFrame for {parser} ({variant}):\n", df.head())

        # (a) トップレイヤーにcanonicalがあるか の円グラフ
        if variant == "all":
            df['is_multilayer'] = df["num_of_layers"] > 1 # 元々 multilayer かどうかの場合分けで、2 通り plotting

            for multilayer in [True, False]:
                sub = df[df['is_multilayer'] == multilayer]
                can_bp_mainlayer = sub['canonical_bp_in_main_layer'] > 0
                non_can_bp_mainlayer = sub['non_canonical_bp_in_main_layer'] > 0
                counts = pd.Series({
                    'Canonical': can_bp_mainlayer.sum(),
                    'Non-Canonical': non_can_bp_mainlayer.sum()
                })

                plt.figure()
                counts.plot.pie(
                    autopct='%1.1f%%',
                    labels=['Canonical', 'Non-Canonical'],
                    startangle=90,
                    title=f'{parser} {"Multilayer" if multilayer else "Single-layer"} ({variant})'
                )
                plt.ylabel('')
                plt.tight_layout()
                fn = f'pie_{"multilayer" if multilayer else "including_singlelayer"}_{variant}.png'
                plt.savefig(os.path.join(output_dir, parser, fn))
                print(f"Saved pie chart for {parser} ({variant}, {'multilayer' if multilayer else 'single-layer'}) to {output_dir / parser / fn}")
                plt.close()
        else: 
            pass

        # (b) pseudoknot layer 数の頻度棒グラフ
        freq = df['num_of_layers'].value_counts(normalize=True).sort_index()
        plt.figure()
        freq.plot.bar()
        plt.xlabel('Pseudoknot Layer Count')
        plt.xticks([int(x) for x in freq.index], rotation=0)
        plt.ylabel('Frequency (%)')
        plt.title(f'{parser} Pseudoknot Layers ({variant})')
        plt.tight_layout()
        fn = f'bar_pseudoknot_layers_{variant}.png'
        plt.savefig(os.path.join(output_dir, parser, fn))
        print(f"Saved bar graph for {parser} ({variant}) to {output_dir / parser / fn}")
        plt.close()