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
                         description='utility script to aggregate summary'
                                     'stats')
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
    parser.add_argument('--diag',
                        default=False,
                        help='creates diagonal plot when `True`')
    args = parser.parse_args()

    df_meta = pd.read_csv(args.input_list, sep='\t')
    df_stat = pd.concat(pd.read_csv(path, sep='\t')
                        for path in df_meta.path)
    assert len(df_meta) == len(df_stat)
    # extract statistic name
    df_stat.reset_index(inplace=True)

    df_stat.dropna(axis=1, how='all', inplace=True)

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

    def draw_heatmap(*these_args, **kwargs):
        '''note: need "these_args" because "args" is already taken'''
        data = kwargs.pop('data').pivot(index=these_args[1],
                                        columns=these_args[0],
                                        values=these_args[2])
        if args.smooth is not None:
            data.loc[:, :] = gaussian_filter(data, sigma=args.smooth)
        sns.heatmap(data, **kwargs).invert_yaxis()

    for metric, group in df.groupby("metric"):
        # n_iid vs. n_epi plot
        fg = sns.FacetGrid(group, col='$d$', height=4)
        fg.map_dataframe(draw_heatmap, '$n_i$', '$n_e$', "value",
                         square=True, vmin=group["value"].min(),
                         # vmax=df.skewness.max(), cbar=False, cmap='Reds',
                         )
        fg.fig.suptitle(metric)
        plt.savefig(f'{args.outdir}/agg_{metric}.pdf')
        group.to_csv(f'{args.outdir}/agg_{metric}.csv')

        if args.diag:
            # prop epi plot
            group = group.assign(**{"proportion epistatic":
                                 lambda x: x['$n_e$'] /
                                 (x['$n_i$'] + x['$n_e$'])})
            max_iid = group['$n_i$'].max()
            group_diag = group.loc[df['$n_i$'] + group['$n_e$'] == max_iid, :]
            plt.figure()
            sns.lmplot(x='proportion epistatic', y="value",
                       data=group_diag, col='$d$')
            plt.xlim([0, 1])
            plt.savefig(f'{args.outdir}/agg_{metric}.diagonal.pdf')


if __name__ == '__main__':
    main()
