#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from Bio import AlignIO
from io import StringIO  # needed to deal with issue NOTE below
from scipy.stats import describe


class AlignMI():
    '''this class augments a biopython MSA with mutual information'''

    def __init__(self, aln_file: str, format='nexus') -> None:
        '''pass an alignment file aln_file
        biopython will hopefully infer the format from the suffix'''
        # NOTE: issue biopython/biopython#1022 can't handle hyphens in taxon
        #       names for nexus parsing. Here's a stupid workaround that parses
        #       line by line and changes hyphens to underscores
        aln_str = ''
        for line in open(aln_file):
            fields = line.split('   ')
            aln_str += '   '.join([fields[0].replace('-', '_'), *fields[1:]])
        self.aln = AlignIO.read(StringIO(aln_str), format)
        # self.aln = AlignIO.read(aln_file, format)
        # END hack
        aln_len = self.aln.get_alignment_length()
        n_taxa = len(self.aln)
        chars = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        n_chars = len(chars)
        # site-wise character counts
        f_single = np.zeros((aln_len, 4))
        for site in range(aln_len):
            for char in self.aln[:, site]:
                f_single[site, chars[char]] += 1 / n_taxa

        # mutual information for all site pairs, right triangular matrix
        self.mi = np.empty((aln_len, aln_len))
        self.mi[:] = np.nan
        for i in range(aln_len):
            for j in range(aln_len):
                if j <= i:
                    continue
                # character pair counts
                f_pair = np.zeros((4, 4))
                for char_i, char_j in zip(self.aln[:, i],
                                          self.aln[:, j]):
                    f_pair[chars[char_i], chars[char_j]] += 1 / n_taxa
                self.mi[i, j] = sum(f_pair[k, l] * np.log(f_pair[k, l]
                                                          / (f_single[i, k]
                                                             * f_single[j, l]))
                                    if f_pair[k, l] > 0 else 0
                                    for k in range(n_chars)
                                    for l in range(n_chars))

    def plot(self, outfile=None) -> plt.Figure:
        '''plot heatmap of MI, optionally writing to file'''
        fig = plt.figure(figsize=(10, 10))
        sns.heatmap(self.mi, square=True)
        if outfile is not None:
            plt.savefig(outfile)
        return fig

    def write(self, outfile) -> None:
        '''write MI matrix to outfile'''
        np.savetxt(outfile, self.mi)

    def write_summary(self, outfile) -> None:
        '''write summary stats of the MI matrix to file'''
        values = self.mi.flatten()
        values = values[~np.isnan(values)]
        stats = describe(values, nan_policy='omit')
        with open(outfile, 'w') as f:
            print('\t'.join(['min', 'max', 'mean', 'variance', 'skewness',
                             'kurtosis']),
                  file=f)
            print(f'{stats.minmax[0]}\t{stats.minmax[1]}\t{stats.mean}\t'
                  f'{stats.variance}\t{stats.skewness}\t{stats.kurtosis}',
                  file=f)


def main():
    '''usage: python align_mi.py -h'''
    import argparse

    parser = argparse.ArgumentParser(
                         description='mutual information of alignment columns')
    parser.add_argument('aln_file',
                        type=str,
                        help='path to alignment file')
    parser.add_argument('--format',
                        type=str,
                        help='alignment file format (default nexus)',
                        default='nexus')
    parser.add_argument('--plot',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    align_mi = AlignMI(args.aln_file, args.format)
    if args.plot:
        align_mi.plot(f'{args.aln_file}.mi.png')
        align_mi.write(f'{args.aln_file}.mi.txt')
    align_mi.write_summary(f'{args.aln_file}.mi.summary.tsv')


if __name__ == '__main__':
    main()
