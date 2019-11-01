#! /usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import os.path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
                         description='calculate pvalues from pps')
    parser.add_argument('input_list',
                        type=str,
                        help='path to file listing paths '
                             'to pps metric values')
    parser.add_argument('aln_metric',
                        type=str,
                        help='path to file with true metric value')
    return parser.parse_args()


def calc_pvalue(pps, true):
    """Calculate pvalue from pps and true metrics.

    >>> calc_pvalue(pd.DataFrame([[1, 2, 3],[4, 5, 6]], \
                                 columns=['m1', 'm2', 'm3']), \
                                 pd.DataFrame([[2, 2, 2],[2, 2, 2]], \
                                 columns=['m1', 'm2', 'm3']))\
                                 .loc[True].tolist()
    [0.5, 1.0, 1.0]
    """
    true = pd.concat([true]*len(pps), ignore_index=True)
    pvalue = pps.subtract(true).ge(0)  # greater than or equal to pps
    pvalue = pvalue.apply(pd.Series.value_counts).fillna(0).loc[True:] / len(pps)
    return pvalue


def main():
    """usage: python pps_pvalue.py -h."""
    # parse command line arguments
    args = parse_args()

    # pps metrics
    with open(args.input_list, 'r') as f:
        fnames = f.read().splitlines()
    df = pd.concat([pd.read_csv(path, sep='\t')
                   for path in fnames])
    reps = len(df)
    df = df.assign(rep=range(reps)).set_index('rep')

    # 'true' metric
    true = pd.read_csv(args.aln_metric, sep='\t')

    # pvalue
    df = calc_pvalue(df, true)

    df = df.rename(columns={i:j for i,j in
                            zip(df.columns.values,
                            [f'{x}_pvalue' for x in df.columns.values])})
    assert df.isnull().sum(axis=1).sum() == 0
    df.to_csv(f'{os.path.splitext(args.aln_metric)[0]}.pvalue.csv',
              index=False, sep='\t')


if __name__ == '__main__':
    main()
