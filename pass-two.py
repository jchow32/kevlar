#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import khmer
from khmer import calc_expected_collisions
from khmer.khmer_args import memory_setting

parser = argparse.ArgumentParser()
parser.add_argument('-M', '--memory', metavar='MEM', type=memory_setting)
parser.add_argument('-k', '--ksize', type=int, default=31, metavar='K')
parser.add_argument('--filters', nargs='+', metavar='FILE')
parser.add_argument('--outfiles', nargs='+', metavar='FILE')
parser.add_argument('--sample', nargs='+', action='append', metavar='FILE',
                    required=True)
args = parser.parse_args()

noutfiles = len(args.outfiles)
nsamples = len(args.sample)
nfilters = len(args.filters)
if noutfiles != nsamples or noutfiles != nfilters:
    message = 'number of outfiles ({:d})'.format(noutfiles)
    message += ' does not match number of samples ({:d})'.format(nsamples)
    message += ' and/or number of Bloom filters ({:d})'.format(nfilters)
    raise IndexError(message)

tablesize = args.memory / 4
counts = [khmer.Counttable(args.ksize, tablesize, 4) for _ in range(nsamples)]
bloomfilters = [khmer.Nodetable.load(f) for f in args.filters]
print('DEBUG loaded')
samples = zip(args.sample, bloomfilters, counts, args.outfiles)
for infilelist, bloomfilter, counttable, outfile in samples:
    for infile in infilelist:
        print('DEBUG infile', infile)
        parser = khmer.ReadParser(infile)
        for n, read in enumerate(parser, 1):
            if n % 10000 == 0:
                print('DEBUG read', read.sequence, n)
            for k, kmer in enumerate(bloomfilter.get_kmers(read.sequence), 1):
                #if k % 1000 == 0:
                #print('DEBUG kmer', n, k, kmer)
                kmer_samples = sum([1 for ct in counts if ct.get(kmer)])
                if kmer_samples == 1:
                    counttable.add(kmer)

    fpr = calc_expected_collisions(counttable, force=True, max_false_pos=1.0)
    print('Estimated FPR: {:.3f}'.format(fpr))
    counttable.save(outfile)
