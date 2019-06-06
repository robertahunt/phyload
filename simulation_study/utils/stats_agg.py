#! /usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.ndimage.filters import gaussian_filter


def main():
    '''usage: python stats_agg.py -h'''
    import argparse

    parser = argparse.ArgumentParser(
                         description='utility script to aggregate mutual '
                                     'information summary stats')
    parser.add_argument('input_list',
                        type=str,
                        help='path to file listing paths '
                             'with their parameters')
    parser.add_argument('outdir',
                        type=str,
                        help='output directory')
    parser.add_argument('--smooth',
                        type=float,
                        default=None,
                        help='gaussian smooth the heatmap data with this '
                             'bandwidth')
    args = parser.parse_args()

    df_meta = pd.read_csv(args.input_list, sep='\t')
    df_stat = pd.concat(pd.read_csv(path, sep='\t')
                        for path in df_meta.path)
    # extract statistic name
    metric = df_stat.columns.values
    assert len(metric) == 1, f'Expected one column, found {len(metric)}'
    metric = metric[0]
    df_stat.reset_index(inplace=True)

    # note the groupby and mean below will average over replicates
    df = pd.concat((df_meta, df_stat),
                   axis=1).groupby(['d',
                                    'n_iid',
                                    'n_epi']).mean().reset_index()

    # fancier TeX column names for prettier plots
    df.rename(index=str,
              columns={'d': '$d$', 'n_iid': '$n_i$', 'n_epi': '$n_e$'},
              inplace=True)

    def draw_heatmap(*these_args, **kwargs):
        '''note: need "these_args" because "args" is already taken'''
        data = kwargs.pop('data').pivot(index=these_args[1],
                                        columns=these_args[0],
                                        values=these_args[2])
        if args.smooth is not None:
            data.loc[:, :] = gaussian_filter(data, sigma=args.smooth)
        sns.heatmap(data, **kwargs).invert_yaxis()
    fg = sns.FacetGrid(df, col='$d$', height=4)
    fg.map_dataframe(draw_heatmap, '$n_i$', '$n_e$', metric,
                     square=True,
                     # vmin=df[metric].min(),
                     # vmax=df.skewness.max(), cbar=False, cmap='Reds',
                     )
    fg.fig.suptitle(metric)
    plt.savefig(f'{args.outdir}/agg_{metric}.pdf')

    df['proportion epistatic'] = df['$n_e$'] / (df['$n_i$'] + df['$n_e$'])
    max_iid = df['$n_i$'].max()
    df_diag = df.loc[df['$n_i$'] + df['$n_e$'] == max_iid, :]
    plt.figure()
    sns.lmplot(x='proportion epistatic', y=metric, data=df_diag, col='$d$')
    plt.xlim([0, 1])
    plt.savefig(f'{args.outdir}/agg_{metric}.diagonal.pdf')




if __name__ == '__main__':
    main()
