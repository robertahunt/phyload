#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from Bio import AlignIO
from io import StringIO # needed to deal with issue NOTE below

class AlignMI():
    '''
    this class augments a biopython MSA with mutual information, uh, information
    '''
    def __init__(self, aln_file: str, format='nexus'):
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
        ### END hack
        variant_sites = [i for i in range(self.aln.get_alignment_length()) if len(set(self.aln[:, i])) > 1]
        n_variant = len(variant_sites)
        n_taxa = len(self.aln)
        chars = {'A':0, 'C':1, 'G':2, 'T':3}
        n_chars = len(chars)
        # site-wise character counts
        f_single = np.zeros((n_variant, 4))
        for i, site in enumerate(variant_sites):
            for char in self.aln[:, site]:
                f_single[i, chars[char]] += 1 / n_taxa

        # mutual information for all site pairs, right triangular matrix
        self.mi = np.empty((n_variant, n_variant))
        self.mi[:] = np.nan
        for i, sitei in enumerate(variant_sites):
            for j, sitej in enumerate(variant_sites):
                if sitej <= sitei: continue
                # character pair counts
                f_pair = np.zeros((4, 4))
                for char_i, char_j in zip(self.aln[:, sitei], self.aln[:, sitej]):
                    f_pair[chars[char_i], chars[char_j]] += 1 / n_taxa
                self.mi[i, j] = sum(f_pair[k, l] * np.log(f_pair[k, l]
                                                          / (f_single[i, k] * f_single[j, l]))
                                    if f_pair[k, l] > 0 else 0
                                    for k in range(n_chars)
                                    for l in range(n_chars))


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
    args = parser.parse_args()

    AlignMI(args.aln_file, args.format)


if __name__  ==  '__main__':
    main()
