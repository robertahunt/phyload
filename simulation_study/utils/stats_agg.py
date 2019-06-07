#! /usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


def main():
    '''usage: python stats_agg.py -h'''
    import argparse

    parser = argparse.ArgumentParser(
                         description='utility script to aggregate summary'
                                     'stats')
    parser.add_argument('input_list',
                        type=str,
                        help='path to file listing paths '
                             'with their parameters')
    parser.add_argument('outdir',
                        type=str,
                        help='output directory')
    args = parser.parse_args()

    df_meta = pd.read_csv(args.input_list, sep='\t')
    df_stat = pd.concat(pd.read_csv(path, sep='\t')
                        for path in df_meta.path)
    # extract statistic name
    df_stat.reset_index(inplace=True)

    # note the groupby and mean below will average over replicates
    df = pd.concat((df_meta, df_stat),
                   axis=1).groupby(['d',
                                    'n_iid',
                                    'n_epi']).mean().reset_index()
    df = pd.melt(df, id_vars=['d', 'n_iid', 'n_epi', 'index'],
                 var_name="metric", value_name="value")

    # fancier TeX column names for prettier plots
    df.rename(index=str,
              columns={'d': '$d$', 'n_iid': '$n_i$', 'n_epi': '$n_e$'},
              inplace=True)

    def draw_heatmap(*args, **kwargs):
        data = kwargs.pop('data').pivot(index=args[1], columns=args[0],
                                        values=args[2])
        sns.heatmap(data, **kwargs).invert_yaxis()
    for metric, group in df.groupby("metric"):
        print(group["value"].min())
        fg = sns.FacetGrid(group, col='$d$', height=4)
        fg.map_dataframe(draw_heatmap, '$n_i$', '$n_e$', "value",
                         square=True, vmin=group["value"].min(),
                         # vmax=df.skewness.max(), cbar=False, cmap='Reds',
                         )
        fg.fig.suptitle(metric)
        plt.savefig(f'{args.outdir}/agg_{metric}.pdf')


if __name__ == '__main__':
    main()
