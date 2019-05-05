#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from glob import glob

def main():
    '''usage: python mi_agg.py -h'''
    import argparse

    parser = argparse.ArgumentParser(
                         description='utility script to aggregate mutual '
                                     'information summary stats')
    parser.add_argument('nest_dir',
                        type=str,
                        help='simulation nest directory, hierarchically '
                             'ordered by d, n_iid, n_epi, replicate')
    args = parser.parse_args()

    paths = glob(f'{args.nest_dir}/*/*/*/*/aln.nex.mi.summary.txt')
    d, n_iid, n_epi = zip(*[path.split('/')[1:4] for path in paths])

    df_meta = pd.DataFrame({'$d$': d, '$n_i$': n_iid,
                            '$n_e$': n_epi})
    df_mi = pd.concat(pd.read_csv(path, sep='\t')
                      for path in paths).reset_index()
    # note the groupby and mean below will average over replicates
    df = pd.concat((df_meta, df_mi),
                   axis=1).groupby(['$d$',
                                    '$n_i$',
                                    '$n_e$']).mean().reset_index()

    def draw_heatmap(*args, **kwargs):
        data = kwargs.pop('data').pivot(index=args[1], columns=args[0],
                                        values=args[2])
        sns.heatmap(data, **kwargs).invert_yaxis()
    fg = sns.FacetGrid(df, col='$d$', height=4)
    fg.map_dataframe(draw_heatmap, '$n_i$', '$n_e$', 'skewness',
                     square=True, vmin=df.skewness.min(),
                     # vmax=df.skewness.max(), cbar=False, cmap='Reds',
                     )
    plt.savefig(f'{args.nest_dir}/mi_agg.skewness.pdf')


if __name__ == '__main__':
    main()
