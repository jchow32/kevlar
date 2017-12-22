#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from math import log10
import re
import sys
import khmer
import kevlar
from kevlar.vcf import Variant
import scipy


def filter_refr(akmers, rkmers, refr):
    """
    Remove k-mers that are present in the reference.

    The `akmers` and `rkmers` variables are lists storing paired k-mers that
    contain an alternate and reference allele, respectively. Each pair of
    k-mers is kept only if the alternate allele k-mer is not present
    (elsewhere) in the genome.
    """
    newa, newr = list(), list()
    for a, r in zip(akmers, rkmers):
        if refr.get(a) == 0:
            newa.append(a)
            newr.append(r)
    return newa, newr


def get_abundances(kmers, casecounts, controlcounts):
    """
    Create a nested list of k-mer abundances.

    Given a list of k-mers, store the k-mer abundances using one sub-list per
    sample, as shown below.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from case/proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]
    """
    abundances = list()
    for _ in range(len(controlcounts) + 1):
        abundances.append(list())

    for kmer in kmers:
        a = casecounts.get(kmer)
        abundances[0].append(a)
        for i in range(len(controlcounts)):
            a = controlcounts[i].get(kmer)
            abundances[i+1].append()

    return abunds


def set_error_rates(error, nsamples):
    if isinstance(error, float):
        errors = [error] * nsamples
    elif isinstance(error, list):
        assert len(error) == nsamples
        for e in error:
            assert isinstance(e, float)
        errors = error
    else:
        message = 'variable {} doesn\'t quack like a float'.format(error)
        message += ' or a list of floats'
        raise ValueError(message)
    return errors


def abund_log_prob(genotype, abundance, mean, sd, error):
    """
    Calculate conditional k-mer abundance probability

    Compute the log (base 10) probability of the given k-mer abundance
    conditioned on the given genotype (copy number). The `genotype` variable
    represents the number of assumed allele copies and is one of {0, 1, 2}
    (corresponding to genotypes {0/0, 0/1, and 1/1}). The `mean` and `sd`
    variables describe a normal distribution of observed abundances of k-mers
    with copy number 2. The `error` parameter is the error rate.
    """
    if genotype == 0:
        return abundance * log10(error)
    if genotype == 1:
        p = scipy.stats.norm.cdf(abund, mean / 2, sd / 2)
        return log10(p)
    if genotype == 2:
        p = scipy.stats.norm.cdf(abund, mean, sd)
        return log10(p)


def likelihood_denovo(altabunds, refrabunds, mean=30.0, sd=8.0, error=0.01):
    """
    Compute the likelihood that a variant is de novo.

    The `altabunds` and `refrabunds` variables should store nested list of
    k-mer abundances, structured as shown below.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]

    The `error` variable can be a single float or a list of floats. If it is a
    single float, the same error rate will be applied to all control samples.
    If it is a list of floats, it must contain 1 value per control sample.
    """
    errors = set_error_rates(error, nsamples=len(refrabunds) + 1)
    logsum = 0.0
    for alt, refr in zip(altabunds[0], refrabunds[0]):
        logsum += abund_log_prob(1, alt, mean=mean, sd=sd)
        logsum += abund_log_prob(1, refr, mean=mean, sd=sd)
    for i in range(len(controlcounts)):
        for a, r, e in zip(altabunds[i+1], refrabunds[i+1], errors):
            logsum += abund_log_prob(0, a, error=e)
            logsum += abund_log_prob(2, r, mean=mean, sd=sd)
    return logsum


def likelihood_false(altabunds, error=0.01):
    """
    Compute the likelihood that a variant is false.

    The `altabunds` variable should store nested list of k-mer abundances,
    structured as shown below.

    abunds = [
        [15, 14, 13, 16, 14, 15, 14, 14],  # k-mer abundances from proband
        [0, 0, 1, 0, 2, 10, 0, 0],  # k-mer abundances from parent/control 1
        [0, 1, 1, 0, 1, 0, 2, 0],  # k-mer abundances from parent/control 2
    ]

    The `error` variable can be a single float or a list of floats. If it is a
    single float, the same error rate will be applied to all case/proband and
    control/parental samples. If it is a list of floats, it must contain 1
    value per sample.
    """
    errors = set_error_rates(error, nsamples=len(altabunds))
    logsum = 0.0
    for abundlist, e in zip(altabunds, errors):
        for abund in abundlist:
            logsum += abund_log_prob(0, abund, error=e)
    return logsum


def likelihood_inherited(altabunds, mean=30.0, sd=8.0, error=0.01):
    """
    Compute the likelihood that a variant is inherited.

    There are 15 valid inheritance scenarios, 11 of which (shown below) result
    in the proband carrying at least one copy of the alternate allele. Sum
    probabilities across all these different scenarios to derive an aggregate
    likelihood that the variant is inherited.

    The other likelihood calculations are implemented to handle an arbitrary
    number of controls, but this can only handle trios.
    """
    scenarios = [
        (1, 0, 1), (1, 0, 2),
        (1, 1, 0), (1, 1, 1), (1, 1, 2),
        (1, 2, 0), (1, 2, 1),
        (2, 1, 1), (2, 1, 2),
        (2, 2, 1), (2, 2, 2),
    ]
    errors = set_error_rates(error, nsamples=3)
    logsum = 0.0
    abundances = zip(altabunds[0], altabunds[1], altabunds[2])
    for a_c, a_m, a_f in abundances:
        maxval = None
        for g_c, g_m, g_f in scenarios:
            testsum = calc_log_prob(g_c, a_c, mean, sd, errors[0]) + \
                      calc_log_prob(g_m, a_m, mean, sd, errors[1]) + \
                      calc_log_prob(g_f, a_f, mean, sd, errors[2]) + \
                      log10(1.0 / 15.0)
            if maxval is None or testsum > maxval:
                maxval = testsum
        logsum += maxval
    return log10(15.0 / 11.0) * logsum  # 1 / (11/15)


def compute_likelihoods(variant, casecounts, controlcounts, refr=None,
                        mean=30.0, sd=8.0, error=0.01, caselabel=None,
                        ctrllabels=None):
    window = variant.info['VW']
    refrwindow = variant.info['VW']
    if window is None or refrwindow is None:
        variant.info['DN'] = '-inf'
        return

    altkmers = casecounts.get_kmers(window)
    refrkmers = casecounts.get_kmers(refrwindow)
    if refr:
        altkmers, refrkmers = filter_refr(altkmers, refrkmers, refr)

    altabunds = get_abundances(altkmers, casecounts, controlcounts)
    refrabunds = get_abundances(refrkmers, casecounts, controlcounts)
    samplelabels = [caselabel] + ctrllabels
    for label, abundlist in zip(samplelabels, altabunds):
        abundstr = ','.join([str(a) for a in abundlist])
        variant.format(label, 'AA', abundstr)

    variant.info['FP'] = likelihood_false(altabunds, error=error)
    variant.info['DN'] = likelihood_denovo(altabunds, refrabunds, mean=mean,
                                           sd=sd, error=error)
    return


def local_to_global(localcoord, subseqid):
    match = re.search('(\S+)_(\d+)-(\d+)', subseqid)
    assert match, 'unable to parse subseqid {:s}'.format(subseqid)
    seqid = match.group(1)
    globaloffset = int(match.group(2))
    globalcoord = globaloffset + localcoord
    return seqid, globalcoord


def call_snv(target, query, offset, length, ksize):
    t = target.sequence[offset:offset+length]
    q = query.sequence[:length]
    diffs = [(i, t[i], q[i]) for i in range(length) if t[i] != q[i]]
    if len(diffs) == 0:
        seqid, globalcoord = local_to_global(offset, target.name)
        nocall = Variant(seqid, globalcoord, '.', '.', NC='perfectmatch', CS=q)
        return [nocall]

    snvs = list()
    for diff in diffs:
        minpos = max(diff[0] - ksize + 1, 0)
        maxpos = min(diff[0] + ksize, length)
        window = q[minpos:maxpos]
        refrwindow = t[minpos:maxpos]

        numoverlappingkmers = len(window) - ksize + 1
        kmers = [window[i:i+ksize] for i in range(numoverlappingkmers)]
        refr = diff[1].upper()
        alt = diff[2].upper()
        localcoord = offset + diff[0]
        seqid, globalcoord = local_to_global(localcoord, target.name)
        snv = Variant(seqid, globalcoord, refr, alt, VW=window, RW=refrwindow)
        snvs.append(snv)
    return snvs


def call_deletion(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize - 1
    window = query.sequence[minpos:maxpos]
    minpos += offset
    maxpos += offset + indellength
    refrwindow = target.sequence[minpos:maxpos]

    refr = target.sequence[offset+leftmatch-1:offset+leftmatch+indellength]
    alt = refr[0]
    assert len(refr) == indellength + 1
    localcoord = offset + leftmatch
    seqid, globalcoord = local_to_global(localcoord, target.name)
    var = Variant(seqid, globalcoord - 1, refr, alt, VW=window, RW=refrwindow)
    return [var]


def call_insertion(target, query, offset, ksize, leftmatch, indellength):
    minpos = leftmatch - ksize + 1
    maxpos = leftmatch + ksize + indellength - 1
    window = query.sequence[minpos:maxpos]
    minpos += offset
    maxpos += offset - indellength
    refrwindow = target.sequence[minpos:maxpos]

    alt = query.sequence[leftmatch-1:leftmatch+indellength]
    refr = alt[0]
    insertion = alt[1:]
    assert len(insertion) == indellength
    localcoord = offset + leftmatch
    seqid, globalcoord = local_to_global(localcoord, target.name)
    var = Variant(seqid, globalcoord - 1, refr, alt, VW=window, RW=refrwindow)
    return [var]


def make_call(target, query, cigar, ksize):
    snvmatch = re.search('^(\d+)D(\d+)M(\d+)D$', cigar)
    snvmatch2 = re.search('^(\d+)D(\d+)M(\d+)D(\d+)M$', cigar)
    if snvmatch:
        offset = int(snvmatch.group(1))
        length = int(snvmatch.group(2))
        return call_snv(target, query, offset, length, ksize)
    elif snvmatch2 and int(snvmatch2.group(4)) <= 5:
        offset = int(snvmatch2.group(1))
        length = int(snvmatch2.group(2))
        return call_snv(target, query, offset, length, ksize)

    indelmatch = re.search('^(\d+)D(\d+)M(\d+)([ID])(\d+)M(\d+)D$', cigar)
    indelmatch2 = re.search('^(\d+)D(\d+)M(\d+)([ID])(\d+)M(\d+)D(\d+)M$',
                            cigar)
    if indelmatch:
        offset = int(indelmatch.group(1))
        leftmatch = int(indelmatch.group(2))
        indellength = int(indelmatch.group(3))
        indeltype = indelmatch.group(4)
        callfunc = call_deletion if indeltype == 'D' else call_insertion
        return callfunc(target, query, offset, ksize, leftmatch, indellength)
    elif indelmatch2 and int(indelmatch2.group(7)) <= 5:
        offset = int(indelmatch2.group(1))
        leftmatch = int(indelmatch2.group(2))
        indellength = int(indelmatch2.group(3))
        indeltype = indelmatch2.group(4)
        callfunc = call_deletion if indeltype == 'D' else call_insertion
        return callfunc(target, query, offset, ksize, leftmatch, indellength)

    seqid, globalcoord = local_to_global(0, target.name)
    nocall = Variant(seqid, globalcoord, '.', '.', NC='inscrutablecigar',
                     CS=query.sequence, CIGAR=cigar)
    return [nocall]


def align_both_strands(targetseq, queryseq, match=1, mismatch=2, gapopen=5,
                       gapextend=0):
    cigar1, score1 = kevlar.align(targetseq, queryseq, match, mismatch,
                                  gapopen, gapextend)
    cigar2, score2 = kevlar.align(targetseq, kevlar.revcom(queryseq), match,
                                  mismatch, gapopen, gapextend)

    if score2 > score1:
        cigar = cigar2
        score = score2
        strand = -1
    else:
        cigar = cigar1
        score = score1
        strand = 1
    return cigar, score, strand


def call(targetlist, querylist, match=1, mismatch=2, gapopen=5, gapextend=0,
         ksize=31, casecounts=None, controlcounts=None, refr=None, mu=30.0,
         sigma=8.0, epsilon=0.01, caselabel=None, ctrllabels=None):
    """
    Wrap the `kevlar call` procedure as a generator function.

    The `targetlist` and `querylist` variables should be iterables containig
    sequences to be aligned. These lists will be sorted and iterated over
    multiple times, so they cannot be used in a streaming fashion.

    The `match, `mismatch`, `gapopen`, and `gapextend` parameters are for
    alignment scoring and should all be non-negative integers.

    If `casecounts` and `controlcounts` are provided, the likelihood of each
    variant will be computed. `casecounts` should be a Counttable with k-mer
    abundances for the proband sample, and `controlcounts` should be a list of
    Counttables, one per parent/control sample. If `refr` is provided (a sketch
    with k-mers from the reference genome), any k-mers containing the alternate
    allele that occur elsewhere in the genome will be discarded.

    The `mu` and `sigma` parameters refer to the mean and standard deviation of
    the observed k-mer abundance distributions, and each can be specified using
    a single float or a list of floats (1 per sample). The `epsilon` parameter
    is an error rate, and can also be specified using a float or a list of
    floats.
    """
    dolike = casecounts is not None and controlcounts is not None
    varcalls = list()
    for query in sorted(querylist, reverse=True, key=len):
        bestcigar = None
        bestscore = None
        besttarget = None
        bestorientation = None
        for target in sorted(targetlist, key=lambda record: record.name):
            cigar, score, strand = align_both_strands(
                target.sequence, query.sequence, match, mismatch, gapopen,
                gapextend
            )
            if bestscore is None or score > bestscore:
                bestscore = score
                bestcigar = cigar
                besttarget = target
                bestorientation = strand

        if bestorientation == -1:
            query.sequence = kevlar.revcom(query.sequence)
        for variant in make_call(besttarget, query, bestcigar, ksize):
            if dolike:
                compute_likelihoods(
                    variant, casecounts, controlcounts, refr=refr, mean=mu,
                    sd=sigma, error=epsilon, caselabel=caselabel,
                    ctrllables=ctrllabels
                )
        varcalls.append(variant)

    if dolike:
        varcalls.sort(key=lambda v: float(v.info['DN']), reverse=True)
    for v in varcalls:
        yield v


def main(args):
    outstream = kevlar.open(args.out, 'w')
    writer = kevlar.vcf.VCFWriter(source='kevlar::call')

    qinstream = kevlar.parse_augmented_fastx(kevlar.open(args.queryseq, 'r'))
    queryseqs = [record for record in qinstream]
    targetseqs = [record for record in khmer.ReadParser(args.targetseq)]

    caselabel = None
    ctrllabels = None
    if args.like_filter:
        refr = khmer.Nodetable.load(args.refr)
        casecounts = khmer.Counttable.load(args.case)
        ctrlcounts = [khmer.Counttable.load(c) for c in args.control]
        caselabel = args.case_label if args.case_label else 'Case'
        writer.register_sample(caselabel)
        if args.ctrl_labels:
            ctrllabels = args.ctrl_labels.split(',')
        else:
            numcontrols = len(args.control)
            ctrllabels = ['Control' + str(i+1) for i in range(numcontrols)]
        for label in ctrllabels:
            writer.register_sample(label)

    caller = call(
        targetseqs, queryseqs, args.match, args.mismatch, args.open,
        args.extend, args.ksize, args.case, args.ctrl, args.refr, args.mu,
        args.sigma, args.epsilon, caselabel, ctrllabels
    )

    for varcall in caller:
        writer.write(varcall, fh=outstream)
