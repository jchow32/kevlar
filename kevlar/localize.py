#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from subprocess import Popen, PIPE
from tempfile import TemporaryFile
import sys

import kevlar
import khmer
import pysam


class KevlarBWAError(RuntimeError):
    """Raised if the delegated BWA call fails for any reason."""
    pass


class KevlarNoReferenceMatchesError(ValueError):
    """Raised if contigs have no k-mer matches against the reference."""
    pass


class KevlarVariantLocalizationError(ValueError):
    """Raised if k-mers match to mutliple locations in the reference."""
    pass


class KevlarRefrSeqNotFound(ValueError):
    """Raised if the reference sequence cannot be found."""
    pass


class IntervalSet(object):

    def __init__(self):
        self._starts = dict()
        self._ends = dict()

    def __len__(self):
        return len(self._starts)

    def add(self, seqid, pos, ksize):
        startpos = max(pos, 0)
        endpos = pos + ksize
        if seqid not in self._starts:
            self._starts[seqid] = startpos
            self._ends[seqid] = endpos
        else:
            self._starts[seqid] = min(startpos, self._starts[seqid])
            self._ends[seqid] = max(endpos, self._ends[seqid])

    def get(self, seqid):
        if seqid not in self._starts:
            return None, None
        return self._starts[seqid], self._ends[seqid]

    @property
    def seqids(self):
        return set(list(self._starts.keys()))


def get_unique_kmers(infile, ksize=31):
    """
    Grab all unique k-mers from the specified sequence file.

    The input file is expected to contain contigs in augmented FASTA format.
    The absence of annotated k-mers is not problematic, but the contigs should
    be on a single line.
    """
    ct = khmer.Counttable(ksize, 1, 1)
    kmers = set()
    instream = kevlar.open(infile, 'r')
    for record in kevlar.parse_augmented_fastx(instream):
        for kmer in ct.get_kmers(record.sequence):
            minkmer = kevlar.revcommin(kmer)
            if minkmer not in kmers:
                kmers.add(minkmer)
                yield kmer


def unique_kmer_string(infile, ksize=31):
    """Convert contigs to k-mer Fasta for BWA input."""
    output = ''
    for n, kmer in enumerate(get_unique_kmers(infile, ksize)):
        output += '>kmer{:d}\n{:s}\n'.format(n, kmer)
    return output


def get_exact_matches(infile, bwaindexfile, ksize=31):
    """
    Compute a list of exact k-mer matches using BWA MEM.

    Input should be a Fasta file containing contigs generated by
    `kevlar assemble`. This function decomposes the contigs into their
    constituent k-mers and searches for exact matches in the reference using
    `bwa mem`. This function is a generator, and yields tuples of
    (seqid, startpos) for each exact match found.
    """
    kmers = unique_kmer_string(infile, ksize)
    cmd = 'bwa mem -k {k} -T {k} {idx} -'.format(k=ksize, idx=bwaindexfile)
    cmdargs = cmd.split(' ')
    with TemporaryFile() as samfile:
        bwaproc = Popen(cmdargs, stdin=PIPE, stdout=samfile, stderr=PIPE,
                        universal_newlines=True)
        stdout, stderr = bwaproc.communicate(input=kmers)
        if bwaproc.returncode != 0:  # pragma: no cover
            print(stderr, file=sys.stderr)
            raise KevlarBWAError('problem running BWA')
        samfile.seek(0)
        sam = pysam.AlignmentFile(samfile, 'r')
        for record in sam:
            if record.is_unmapped:
                continue
            seqid = sam.get_reference_name(record.reference_id)
            yield seqid, record.pos


def extract_regions(refr, seedmatches, maxspan=1000, delta=50):
    """
    Extract the specified genomic region from the provided file object.

    The start and end parameters define a 0-based half-open genomic interval.
    Bounds checking must be performed on the end parameter.
    """
    observed_seqids = set()
    for defline, sequence in kevlar.seqio.parse_fasta(refr):
        seqid = defline[1:].split()[0]
        observed_seqids.add(seqid)

        start, end = seedmatches.get(seqid)
        if start is None:
            continue

        newstart = max(start - delta, 0)
        newend = min(end + delta, len(sequence))
        span = newend - newstart
        if span > maxspan:
            message = 'variant spans {:d} bp (max {:d})'.format(span, maxspan)
            raise KevlarVariantLocalizationError(message)

        subseqid = '{}_{}-{}'.format(seqid, newstart, newend)
        subseq = sequence[newstart:newend]
        yield subseqid, subseq

    missing = [s for s in seedmatches.seqids if s not in observed_seqids]
    if len(missing) > 0:
        raise KevlarRefrSeqNotFound(','.join(missing))


def main(args):
    instream = kevlar.open(args.refr, 'r')
    output = kevlar.open(args.out, 'w')

    seedmatches = IntervalSet()
    for seqid, pos in get_exact_matches(args.contigs, args.refr, args.ksize):
        seedmatches.add(seqid, pos, args.ksize)
    if len(seedmatches) == 0:
        raise KevlarNoReferenceMatchesError()

    for subseqid, subseq in extract_regions(instream, seedmatches):
        print('>', subseqid, '\n', subseq, sep='', file=output)
