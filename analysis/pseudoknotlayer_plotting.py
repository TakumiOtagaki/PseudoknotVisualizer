import os
import json
import pandas as pd
import matplotlib.pyplot as plt

# 設定
parsers = ['rnaview', 'dssr']
base_dir = 'analysis'
output_dir = 'analysis/graphs'

# 出力先ディレクトリ作成
for parser in parsers:
    os.makedirs(os.path.join(output_dir, parser), exist_ok=True)

for parser in parsers:
    for variant in ['canonical_only', 'all']:
        # JSON読み込み
        path = os.path.join(base_dir, f'pseudoknot_analysis_{parser}.{variant}.json')
        with open(path) as f:
            data = json.load(f)
        df = pd.DataFrame(data)

        # (a) トップレイヤーにcanonicalがあるか の円グラフ
        df['is_multilayer'] = df['pseudoknot_layer_count'] > 1
        def top_has_canonical(layers):
            top = layers[-1]
            return any(bp['is_canonical'] for bp in top['basepair_details'])
        df['top_canonical'] = df['layers'].apply(top_has_canonical)

        for multilayer in [True, False]:
            sub = df[df['is_multilayer'] == multilayer]
            counts = sub['top_canonical'].value_counts(normalize=True)
            plt.figure()
            counts.plot.pie(
                autopct='%1.1f%%',
                labels=['No Canonical', 'Has Canonical'],
                startangle=90,
                title=f'{parser} {"Multilayer" if multilayer else "Single-layer"} ({variant})'
            )
            plt.ylabel('')
            plt.tight_layout()
            fn = f'pie_{"multilayer" if multilayer else "singlelayer"}_{variant}.png'
            plt.savefig(os.path.join(output_dir, parser, fn))
            plt.close()

        # (b) pseudoknot layer 数の頻度棒グラフ
        freq = df['pseudoknot_layer_count'].value_counts(normalize=True).sort_index()
        plt.figure()
        freq.plot.bar()
        plt.xlabel('Pseudoknot Layer Count')
        plt.ylabel('Frequency (%)')
        plt.title(f'{parser} Pseudoknot Layers ({variant})')
        plt.tight_layout()
        fn = f'bar_pseudoknot_layers_{variant}.png'
        plt.savefig(os.path.join(output_dir, parser, fn))
        plt.close()