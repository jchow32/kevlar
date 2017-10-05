#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import threading
import khmer
from khmer.khmer_args import memory_setting

parser = argparse.ArgumentParser()
parser.add_argument('-M', '--memory', metavar='MEM', type=memory_setting)
parser.add_argument('-k', '--ksize', type=int, default=31, metavar='K')
parser.add_argument('-t', '--threads', type=int, default=1, metavar='T')
parser.add_argument('outfile')
parser.add_argument('infiles', nargs='+')
args = parser.parse_args()

tablesize = args.memory * khmer._buckets_per_byte['nodegraph'] / 4
nt = khmer.Nodetable(args.ksize, tablesize, 4)
for infile in args.infiles:
    parser = khmer.ReadParser(infile)
    threads = list()
    for _ in range(args.threads):
        t = threading.Thread(target=nt.consume_seqfile, args=(parser, ))
        threads.append(t)
        t.start()
    for t in threads:
        t.join()

fpr = khmer.calc_expected_collisions(nt, force=True, max_false_pos=1.0)
print('Estimated FPR: {:.3f}'.format(fpr))
nt.save(args.outfile)
