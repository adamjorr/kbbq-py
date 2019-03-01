#!/usr/bin/env python3
import pysam
import numpy as np
import sklearn
from sklearn.linear_model import LogisticRegression as LR
from sklearn.isotonic import IsotonicRegression as IR
import os
import os.path
import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import khmer
import scipy.stats
from . import recaltable
import pandas as pd

def tstamp():
    return '[ ' + datetime.datetime.today().isoformat(' ', 'seconds') + ' ]'

def load_positions(posfile):
    d = dict()
    with open(posfile, 'r') as infh:
        for line in infh:
            chrom, pos = line.rstrip().split()
            d.setdefault(chrom, list()).append(int(pos))
    return d

def find_rcorrected_sites(uncorrfile, corrfile):
    print(tstamp(), "Finding Rcorrected sites . . .", file=sys.stderr)
    uncorr_reads = list(pysam.FastxFile(uncorrfile))
    corr_reads = list(pysam.FastxFile(corrfile))
    #verify the sequences are the same and can be accessed by index
    for i in range(len(corr_reads)):
        try:
            assert corr_reads[i].name.startswith(uncorr_reads[i].name)
        except AssertionError:
            print("Corr_set[i]:",corr_reads[i])
            print("Uncorr_Set[i]:", uncorr_reads[i])
            raise

    names = dict()
    seqlen = len(uncorr_reads[0].get_quality_array())
    rawquals = np.zeros([len(uncorr_reads), seqlen], dtype = np.int)
    rcorrected = np.zeros([len(uncorr_reads),seqlen], dtype = np.bool)
    seqs = np.zeros(len(uncorr_reads), dtype = 'U' + str(seqlen))
    rgs = np.zeros(len(uncorr_reads), dtype = 'U7')
    for i in range(len(uncorr_reads)):
        names[uncorr_reads[i].name.split(sep='_')[0]] = i
        #print(uncorr_reads[i].name)
        rgs[i] = uncorr_reads[i].name.split(sep='_')[1].split(':')[-1]
        rawquals[i,:] = uncorr_reads[i].get_quality_array()
        seqs[i] = uncorr_reads[i].sequence
        uncorr_s = np.array(list(uncorr_reads[i].sequence), dtype = np.unicode)
        corr_s = np.array(list(corr_reads[i].sequence), dtype = np.unicode)
        rcorrected[i] = (uncorr_s == corr_s)
    return names, rawquals.copy(), rcorrected.copy(), seqs.copy(), rgs.copy(), seqlen

def train_regression(rawquals, rcorrected, tol = 1e-4):
    print(tstamp(), "Doing Logit Regression", file=sys.stderr)
    lr = LR(tol = tol)
    lr = lr.fit(rawquals.flatten().reshape(-1,1), rcorrected.flatten())
    return lr

def recalibrate(lr, q):
    print(tstamp(), "Recalibrating Quality Scores . . .", file = sys.stderr)
    newprobs = lr.predict_proba(q.flatten().reshape(-1,1))[:,1]
    newq = p_to_q(newprobs)
    newq.reshape(q.shape)
    assert newq.shape == q.shape
    return q.copy()

def find_read_errors(read, ref, variable):
    #use the CIGAR to find errors in the read or sites to skip
    #we will add softclipped bases to a skip array and return the error array and the skip array
    #we don't consider indel errors and just track them to properly navigate the reference
    # here's how gatk does it: https://github.com/broadinstitute/gatk/blob/78df6b2f6573b3cd2807a71ec8950d7dfbc9a65d/src/main/java/org/broadinstitute/hellbender/utils/recalibration/BaseRecalibrationEngine.java#L370
    seq = np.array(list(read.query_sequence), dtype = np.unicode)
    skips = np.zeros(seq.shape, dtype = np.bool)
    cigartuples = read.cigartuples #list of tuples [(operation, length)]
    cigarops, cigarlen = zip(*cigartuples)
    cigarops = np.array(cigarops, dtype = np.int)
    cigarlen = np.array(cigarlen, dtype = np.int)

    #reference length from CIGAR: https://github.com/samtools/htsjdk/blob/942e3d6b4c28a8e97c457dfc89625bb403bdf83c/src/main/java/htsjdk/samtools/Cigar.java#L76
    #sum lengths of MDN=X
    #reflen = np.sum(cigarlen[np.any([cigarops == 0, cigarops == 2, cigarops == 3, cigarops == 7, cigarops == 8], axis = 0)])
    subset_variable = variable[read.reference_name][read.reference_start : read.reference_end]
    refseq = np.array(list(ref[read.reference_name][read.reference_start : read.reference_end]), dtype = np.unicode)

    readerrors = np.zeros(seq.shape, dtype = np.bool)
    readidx = 0
    refidx = 0
    for op, l in cigartuples:
        if op == 0 or op == 7 or op == 8:
            #match
            readerrors[readidx : readidx + l] = (refseq[refidx : refidx + l] != seq[readidx : readidx + l])
            skips[readidx : readidx + l] = subset_variable[refidx : refidx + l]
            readidx = readidx + l
            refidx = refidx + l
        elif op == 1:
            #insertion in read
            skips[readidx : readidx + l] = True
            readidx = readidx + l
        elif op == 2 or op == 3:
            #deletion in read or N op
            # N is for introns in mRNA
            refidx = refidx + l
        elif op == 4:
            #soft clip, consumes query not ref
            skips[readidx : readidx + l] = True
            readidx = readidx + l
        elif op == 5 or op == 6:
            #hard clip or pad, do nothing
            continue
        else:
            #unrecognized
            raise ValueError("Uncrecognized Cigar Operation " + str(op) + " In Read\n" + str(read))
    return readerrors, skips

def find_errors(bamfilename, fastafilename, var_pos, names, seqlen):
    #this function may be better optimized that using the pileup
    #since we have to jump around a lot when using the pileup method
    #need to return gatkcalibratedquals, erroneous, skips
    print(tstamp(), "Finding Errors...", file = sys.stderr)
    #rawquals = np.zeros([len(names), seqlen], dtype = np.int)
    gatkcalibratedquals = np.zeros([len(names), seqlen], dtype = np.int)
    erroneous = np.zeros([len(names), seqlen], dtype = np.bool)
    #seqs = np.zeros([len(names), seqlen], dtype = np.unicode)
    skips = np.zeros([len(names),seqlen], dtype = np.bool)
    fasta = pysam.FastaFile(fastafilename)
    ref = {chrom : np.array(list(fasta.fetch(reference = chrom)), dtype = np.unicode) for chrom in fasta.references}
    varsites = {chrom : np.array(var_pos[chrom], dtype = np.int) for chrom in var_pos.keys()}
    fullskips = {chrom : np.zeros(len(ref[chrom]), dtype = np.bool) for chrom in ref.keys()}
    for chrom in fullskips.keys():
        variable_positions = varsites[chrom]
        fullskips[chrom][variable_positions] = True

    bam = pysam.AlignmentFile(bamfilename, 'r')
    readcounter = 0
    for read in bam:
        suffix = ("/2" if read.is_read2 else "/1")
        readidx = names.get(read.query_name + suffix)
        if readidx is None:
            continue

        gatkcalibratedquals[readidx,:] = np.array(read.query_qualities, dtype = np.int)
        e, s = find_read_errors(read, ref, fullskips)
        erroneous[readidx,:] = e
        skips[readidx,:] = s

        readcounter = readcounter + 1

        if read.is_reverse:
            gatkcalibratedquals[readidx,:] = np.flip(gatkcalibratedquals[readidx,:])
            erroneous[readidx,:] = np.flip(erroneous[readidx,:])
            skips[readidx,:] = np.flip(skips[readidx,:])

    try:
        assert readcounter == len(names)
    except AssertionError:
        print("readcounter",readcounter)
        print("len(names)",len(names))
        print("num reads missing:", np.sum(np.all(gatkcalibratedquals == 0, axis = 1)))
        raise
    return gatkcalibratedquals, erroneous, skips

class RescaledNormal:
    oldset = np.seterr(all = 'raise')
    maxscore = 42
    possible_diffs = np.arange(maxscore+1, dtype = np.int_)
    prior_dist = np.zeros(possible_diffs.shape[0], dtype = np.longdouble)
    for i in range(possible_diffs.shape[0]):
        try:
            prior_dist[i] = np.log(.9 * np.exp(-((possible_diffs[i]/.5)**2)/2))
        except FloatingPointError:
            prior_dist[i] = np.NINF
    np.seterr(**oldset)

    def prior(self, difference):
        return prior_dist[difference]

def gatk_delta_q(prior_q, numerrs, numtotal, maxscore = 42):
    assert prior_q.shape == numerrs.shape == numtotal.shape
    possible_q = np.arange(maxscore, dtype = np.int)
    diff = np.absolute(np.subtract.outer(possible_q, prior_q).astype(np.int64))
    #1st dim is possible qs
    prior = RescaledNormal.prior_dist[diff]
    broadcast_errs = np.broadcast_to(numerrs, possible_q.shape + numerrs.shape).copy()
    broadcast_tot = np.broadcast_to(numtotal, possible_q.shape + numtotal.shape).copy()
    p = q_to_p(possible_q).astype(np.float)
    while len(p.shape) < len(broadcast_tot.shape):
        p = np.expand_dims(p, -1)
    broadcast_p = np.broadcast_to( p, broadcast_tot.shape ).copy()
    loglike = scipy.stats.binom.logpmf(broadcast_errs+1, broadcast_tot+2, broadcast_p)
    #loglike should now have same dims as prior
    assert loglike.shape == prior.shape
    posterior = prior + loglike
    posterior_q = np.argmax(posterior, axis = 0)
    try:
        assert posterior_q.shape == prior_q.shape
    except AssertionError:
        print("Posterior",posterior)
        print("Posterior.shape",posterior.shape)
        print("Posterior q:",posterior_q)
        print("Posterior q.shape",posterior_q.shape)
        print("Prior q:",prior_q)
        print("Prior q.shape:", prior_q.shape)
        raise
    return posterior_q - prior_q

#this passes bam test
def table_to_vectors(tablefile, rg_order, dinuc_order, seqlen, maxscore = 42):
    #the recal table uses the PU of the read group as the read group entry in the table
    table = recaltable.RecalibrationReport.from_file(tablefile)
    rgtable = table.tables[2].data.reindex(rg_order)
    meanq = rgtable['EstimatedQReported'].values.astype(np.float64)
    global_errs = rgtable['Errors'].values.astype(np.int64)
    global_total = rgtable['Observations'].values

    qtable = table.tables[3].data.reindex(pd.MultiIndex.from_product([rg_order, np.arange(maxscore + 1)]))
    q_shape = (len(rg_order), maxscore + 1)
    q_errs = qtable['Errors'].fillna(0, downcast = 'infer').values.reshape(q_shape)
    q_total = qtable['Observations'].fillna(0, downcast = 'infer').values.reshape(q_shape)

    postable = table.tables[4].data.loc[rg_order, np.arange(maxscore + 1), 'Cycle']
    postable = postable.reset_index(level = 'CovariateValue').astype({'CovariateValue' : np.int_}).set_index('CovariateValue', append = True)
    postable = postable.reindex(pd.MultiIndex.from_product(
        [rg_order, np.arange(maxscore + 1), ['Cycle'],
        np.concatenate([np.arange(seqlen)+1, np.flip(-(np.arange(seqlen)+1),axis = 0)])
        ]))
    pos_shape = (len(rg_order), maxscore + 1, 2 * seqlen)
    pos_errs = postable['Errors'].fillna(0, downcast = 'infer').values.reshape(pos_shape)
    pos_total = postable['Observations'].fillna(0, downcast = 'infer').values.reshape(pos_shape)

    dinuctable = table.tables[4].data.reindex(pd.MultiIndex.from_product([rg_order, np.arange(maxscore + 1), ['Context'], dinuc_order]))
    dinuc_shape = (len(rg_order), maxscore + 1, len(dinuc_order))
    dinuc_errs = dinuctable['Errors'].fillna(0, downcast = 'infer').values.reshape(dinuc_shape)
    dinuc_total = dinuctable['Observations'].fillna(0, downcast = 'infer').values.reshape(dinuc_shape)

    return meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

def vectors_to_table(meanq, global_errs, global_total, q_errs, q_total,
    pos_errs, pos_total, dinuc_errs, dinuc_total, rg_order, dinuc_order,
    seqlen, maxscore = 42):
    """
    Turn the set of recalibration vectors into a
    :class:`kbbq.recaltable.RecalibrationReport` object.

    For the recalibration vectors, each dimension corresponds to a covariate.
    The first index is always the read group, and the second (if it exists)
    represents the raw quality score, the final index is either the cycle or
    dinucleotide covariate.

    :param np.array[\:] meanq: Mean q for each read group
    :param np.array[\:] global_errs: Number of errors for each read group
    :param np.array[\:] global_total: Number of observations for each read group
    :param np.array[\:,\:] q_errs: Number of errors for each read group and q
        score subset.
    :param np.array[\:,\:] q_total: Number of observations for each read group
        and q score subset.
    :param np.array[\:,\:,\:] pos_errs: Number of errors for each read group, q,
        and cycle subset.
    :param np.array[\:,\:,\:] pos_total: Number of observations for each read
        group, q, and cycle subset.
    :param np.array[\:,\:,\:] dinuc_errs: Number of errors for each read group, q,
        and dinucleotide subset.
    :param np.array[\:,\:,\:] dinuc_total: Number of observations for each read
        group, q, and dinucleotide subset.
    :param list(str) rg_order: The order of read groups
    :param list(str) dinuc_order: The order of dinucleotides
    :param int seqlen: The length of the sequences
    :param int maxscore: The maximum possible quality score
    :return: the recalibration table
    :rtype: :class:`kbbq.recaltable.RecalibrationReport`

    """
    rgdata = {'ReadGroup' : rg_order,
        'EventType' : 'M',
        'EmpiricalQuality' : gatk_delta_q(meanq, global_errs, global_total) + meanq,
        'EstimatedQReported' : meanq,
        'Observations' : global_total,
        'Errors' : global_errs
        }
    rgtable = pd.DataFrame(data = rgdata, index = 'ReadGroup')

    pass

def table_recalibrate(q, tablefile, rg_order, dinuc_order, seqlen, reversecycle, rgs, dinucleotide, minscore = 6, maxscore = 42):
    meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total = table_to_vectors(tablefile, rg_order, dinuc_order, seqlen, maxscore)
    globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq = get_delta_qs(meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total)

    recal_q = np.array(q, copy = True, dtype = np.int)
    valid_positions = (q >= minscore)
    pos = np.broadcast_to(np.arange(q.shape[1]), (q.shape[0], q.shape[1])).copy()
    np.add.at(pos, reversecycle, 1)
    np.negative.at(pos,reversecycle)
    rgcov = rgs[valid_positions]
    qcov = q[valid_positions]
    poscov = pos[valid_positions]
    dinuccov = dinucleotide[valid_positions]

    recal_q[valid_positions] = (meanq[rgcov] + globaldeltaq[rgcov] + qscoredeltaq[rgcov,qcov] + dinucdeltaq[rgcov, qcov, dinuccov] + positiondeltaq[rgcov, qcov, poscov]).astype(np.int)
    return recal_q


def delta_q_recalibrate(q, rgs, dinucleotide, errors, reversecycle, minscore = 6, maxscore = 42):
    print(tstamp(), "Getting Covariate Arrays . . .", file=sys.stderr)
    meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total = get_covariate_arrays(q, rgs, dinucleotide, errors, reversecycle)
    print(tstamp(), "Finding Delta Q's . . .", file=sys.stderr)
    globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq = get_delta_qs(meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total)
    print(tstamp(), "Recalibrating . . .", file=sys.stderr)
    recal_q = np.array(q, copy = True, dtype = np.int)
    
    #vectorize cycle covariates:
    pos = np.broadcast_to(np.arange(q.shape[1]), (q.shape[0], q.shape[1])).copy()
    np.add.at(pos, reversecycle, 1)
    np.negative.at(pos,reversecycle)

    #validate positions with minscore
    valid_positions = (q >= minscore)
    rgcov = rgs[valid_positions]
    qcov = q[valid_positions]
    poscov = pos[valid_positions]
    dinuccov = dinucleotide[valid_positions]

    recal_q[valid_positions] = (meanq[rgcov] + globaldeltaq[rgcov] + qscoredeltaq[rgcov,qcov] + positiondeltaq[rgcov, qcov, poscov] + dinucdeltaq[rgcov, qcov, dinuccov]).astype(np.int)
    return recal_q.copy()

def get_dinucleotide(seqs, q):
    #[A, T, G, C] -> [A, T, G, C]
    #nucleotides at the beginning of the sequence have an empty string before them
    #we should: ignore any context containing an N, ignore any context at beginning of sequence
    #we also ignore the longest string at the start and end of the sequence that have scores <=2
    print(tstamp(), "Getting dinucleotide context . . .", file=sys.stderr)
    nucs = ['A','T','G','C']
    seqs = seqs.copy() #we may need to alter this
    dinucs = [i + j for i in  nucs for j in nucs]
    dinuc_to_int = dict(zip(dinucs, range(len(dinucs))))
    dinucleotide = generic_dinuc_covariate(seqs.view('U1').reshape((seqs.size, -1)), q, dinuc_to_int)
    return dinucleotide.copy(), dinucs.copy()

def get_covariate_arrays(q, rgs, dinucleotide, errors, reversecycle, maxscore = 42, minscore = 6):
    #input arrays are the same dimensions: (numsequences, seqlen)
    #output arrays are dimension (nrgs), (nrgs, q), (nrgs, q, seqlen), or (nrgs, q, 16)
    # m = np.ma.getmaskarray(q)
    nrgs = np.unique(rgs).shape[0]
    seqlen = q.shape[1]

    meanq = np.zeros(nrgs, dtype = np.int_)
    expected_errs = np.zeros(nrgs, dtype = np.longdouble)
    rg_errs = np.zeros(nrgs, dtype = np.int_)
    rg_total = np.zeros(nrgs, dtype = np.int_)
    q_errs = np.zeros((nrgs, maxscore + 1), dtype = np.int_)
    q_total = np.zeros((nrgs, maxscore + 1), dtype = np.int_)
    pos_errs = np.zeros((nrgs, maxscore + 1, 2 * seqlen), dtype = np.int_)
    pos_total = np.zeros((nrgs, maxscore + 1, 2 * seqlen), dtype = np.int_)
    dinuc_errs = np.zeros((nrgs, maxscore + 1, 16), dtype = np.int_)
    dinuc_total = np.zeros((nrgs, maxscore + 1, 16), dtype = np.int_)
    pos = np.broadcast_to(np.arange(seqlen), (q.shape[0], seqlen)).copy()
    np.add.at(pos, reversecycle, 1)
    np.negative.at(pos,reversecycle)

    #these will be reused a lot; cache them here
    # e = np.logical_and(errors, ~m)
    rge = rgs[errors]
    qe = q[errors]
    # valid = np.logical_and(dinucleotide != -1, ~m)
    valid = (dinucleotide != -1)
    e_and_valid = np.logical_and(errors, valid)

    np.add.at(expected_errs, rgs, q_to_p(q))
    np.add.at(rg_errs, rge, 1)
    np.add.at(rg_total, rgs, 1)
    np.add.at(q_errs, (rge, qe), 1)
    np.add.at(q_total, (rgs, q), 1)
    np.add.at(pos_errs, (rge, qe, pos[errors]), 1)
    np.add.at(pos_total, (rgs, q, pos), 1)
    np.add.at(dinuc_errs, (rgs[e_and_valid], q[e_and_valid], dinucleotide[e_and_valid]), 1)
    np.add.at(dinuc_total, (rgs[valid], q[valid], dinucleotide[valid]), 1)

    meanq = p_to_q(expected_errs / rg_total)

    return meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

#this function passes bam test
def get_delta_qs(meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total, maxscore = 42):
    # shapes are:
    #   [rg]
    #   [rg, q]
    #   [rg, q, covariate]
    #seqlen = pos_total.shape[2] / 2
    nrgs = meanq.shape[0]
    rgdeltaq = gatk_delta_q(meanq, rg_errs, rg_total)
    # the qscoredeltaq is 1d, with the index being the quality score
    prior1 = np.broadcast_to((meanq + rgdeltaq)[:,np.newaxis], q_total.shape).copy()
    qscoredeltaq = gatk_delta_q( prior1 , q_errs, q_total)
    ## positiondeltaq is 2d, first dimension is quality score and second is position
    ## dinucdeltaq is 2d, first dimension is quality score and second is nucleotide context
    prior2 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], pos_total.shape).copy()
    positiondeltaq = gatk_delta_q(prior2, pos_errs, pos_total)
    prior3 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], dinuc_total.shape).copy()
    dinucdeltaq = gatk_delta_q(prior3, dinuc_errs, dinuc_total)

    #need to add another value of dinuc, for invalid dinuc
    pad = np.zeros((len(dinucdeltaq.shape),2), dtype = np.int_)
    pad[-1,1] = 1 #add a 0 to the last axis
    dinucdq = np.pad(dinucdeltaq, pad_width = pad, mode = 'constant', constant_values = 0)

    return rgdeltaq.copy(), qscoredeltaq.copy(), positiondeltaq.copy(), dinucdq.copy()

def plot_calibration(data, truth, labels, plotname, plottitle = None, maxscore = 42):
    print(tstamp(), "Making Quality Score Plot . . .", file = sys.stderr)
    assert np.ndim(truth) == 1
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    qualplot, (ax1, ax2) = plt.subplots(2, sharex = True, gridspec_kw = { 'height_ratios' : [2, 1] }, figsize = (7,9))
    #ax1.set_aspect('equal')
    # ax2.set_ylim(top = 3e7)
    qualplot.suptitle(plottitle)
    ax1.plot(np.arange(maxscore + 1), 'k:', label = "Perfect")
    try:
        assert np.all([np.array_equal(data[0].shape, data[i].shape) for i in range(len(data))])
    except AssertionError:
        for i in range(len(data)):
            print(labels[i], 'Shape:', data[i].shape)
        raise
    for i in range(len(data)):
        label = labels[i]
        print(tstamp(), "Plotting %s . . ." % (label), file = sys.stderr)
        assert np.ndim(data[i]) == 1
        est_p = q_to_p(data[i])
        try:
            bscore = sklearn.metrics.brier_score_loss(truth.reshape(-1), est_p.reshape(-1))
        except ValueError:
            print(est_p)
            raise
        numtotal = np.bincount(data[i].reshape(-1), minlength = (maxscore+1))
        numerrs = np.bincount(data[i][truth].reshape(-1), minlength = len(numtotal)) #not masked and error
        p = np.true_divide(numerrs[numtotal != 0],numtotal[numtotal != 0])
        q = p_to_q(p)
        ax1.plot(np.arange(len(numtotal))[numtotal != 0], q, 'o-', alpha = .6, label ="%s, %1.5f" % (label, bscore))
        ax2.plot(np.arange(len(numtotal))[numtotal != 0], numtotal[numtotal != 0], 'o-', alpha = .6, label = "%s, %1.5f" % (label, bscore))
    plt.xlabel("Predicted Quality Score")
    ax1.set_ylabel("Actual Quality Score")
    ax1.legend(loc = "upper left")
    ax2.set_ylabel("Sample Size")
    qualplot.savefig(plotname)

def p_to_q(p, maxscore = 42):
    q = np.zeros(p.shape, dtype = np.int)
    q[p != 0] = (-10.0*np.log10(p[p != 0])).astype(np.int) #avoid divide by 0
    q[p == 0] = maxscore
    q = np.clip(q, 0, maxscore)
    return q.copy()

def q_to_p(q):
    p = np.array(np.power(10.0,-(q / 10.0)), dtype = np.longdouble, copy = True)
    return p

def bam_test(bamfile, tablefile, rg_to_int, rg_order, dinuc_order, seqlen, minscore = 6, maxscore = 42):
    print(tstamp(), "Beginning BAM test . . .", file = sys.stderr)
    meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total = table_to_vectors(tablefile, rg_order, dinuc_order, seqlen, maxscore)
    globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq = get_delta_qs(meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total)
    dinuc_to_int = {d : i for i,d in enumerate(dinuc_order)}
    bam = pysam.AlignmentFile(bamfile,'r')
    for read in bam:
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        recalibrated_quals = recalibrate_bamread(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg_to_int, dinuc_to_int)
        try:
            assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)
        except AssertionError:
            print('GATK calibrated:', gatk_calibrated_quals)
            print('Recalibrated:', recalibrated_quals)
            raise
    print(tstamp(), "BAM test completed successfully.", file = sys.stderr)

## Generic covariate functions

def generic_cycle_covariate(sequencelen, secondinpair = False):
    cycle = np.arange(sequencelen)
    if secondinpair:
        cycle = np.negative(cycle + 1)
    return cycle

def generic_dinuc_covariate(sequences, quals, dinuc_to_int, minscore = 6):
    #this should be refactored to compensate for multiple seqs, multiple quals. it takes too long as is
    #sequences should be a numpy unicode character array, quals should be a numpy array of integer quality scores
    assert sequences.shape == quals.shape
    assert sequences.dtype == np.dtype('U1')
    dinuc = np.char.add(sequences[...,:-1], sequences[...,1:])
    dinuccov = np.zeros(sequences.shape, dtype = np.int)
    dinuccov[...,0] = -1
    is_n = (sequences[...,1:] == 'N')
    follows_n = (sequences[...,:-1] == 'N')
    invalid = np.logical_or(quals[...,1:] < minscore, np.logical_or(is_n, follows_n))
    dinuccov[...,1:][...,invalid] = -1
    vecget = np.vectorize(dinuc_to_int.get)
    dinuccov[...,1:][...,np.logical_not(invalid)] = vecget(dinuc[...,np.logical_not(invalid)])
    return dinuccov

## Recalibrating FASTQ reads

def fastq_cycle_covariates(read, secondinpair = False):
    return generic_cycle_covariate(len(read.sequence), secondinpair)

def fastq_dinuc_covariates(read, dinuc_to_int, minscore = 6):
    quals = np.array(read.get_quality_array(), dtype = np.int)
    return generic_dinuc_covariate(read.sequence, quals, dinuc_to_int, minscore)

def recalibrate_fastq(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg, dinuc_to_int, secondinpair = False, minscore = 6, maxscore = 42):
    qcov = np.array(read.get_quality_array(), dtype = np.int)
    recalibrated_quals = np.array(qcov, copy = True, dtype = np.int)
    valid_positions = (qcov >= minscore)
    cycle = fastq_cycle_covariates(read)[valid_positions]
    dinuccov = fastq_dinuc_covariates(read, dinuc_to_int, minscore)[valid_positions]
    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] + qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuccov] + positiondeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

## Recalibrate reads from a BAM

def bamread_get_fwd_oq(read):
    oq = np.array(list(read.get_tag('OQ')), dtype = np.unicode)
    quals = np.array(oq.view(np.uint32) - 33, dtype = np.uint32)
    if read.is_reverse:
        quals = np.flip(quals)
    return quals

def bamread_cycle_covariates(read):
    cycle = generic_cycle_covariate(len(read.query_sequence), read.is_read2)
    if read.is_reverse:
        cycle = np.flip(cycle)
    return cycle

def bamread_dinuc_covariates(read, dinuc_to_int, complement, minscore = 6):
    #TODO: add stuff to check whether OQ is present,
    # otherwise use read.query_qualities
    seq = read.query_sequence
    oq = np.array(list(read.get_tag('OQ')), dtype = np.unicode_)
    quals = np.array(oq.view(np.uint32) - 33, dtype = np.uint32)
    if read.is_reverse:
        seq = ''.join([complement.get(x,'N') for x in reversed(seq)])
        quals = np.flip(quals)
    dinuccov = generic_dinuc_covariate(seq, quals, dinuc_to_int, minscore)
    if read.is_reverse:
        dinuccov = np.flip(dinuccov)
    return dinuccov

def recalibrate_bamread(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg_to_int, dinuc_to_int, minscore = 6, maxscore = 42):
    complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    
    #TODO: add an assert or logic to ensure OQ is present
    # alternatively, use OQ if present otherwise assume
    # read.query_qualities is the raw score
    oq = np.array(list(read.get_tag('OQ')), dtype = np.unicode_)
    original_quals = np.array(oq.view(np.uint32) - 33, dtype = np.uint32)
    recalibrated_quals = np.array(original_quals, dtype = np.int)
    rg = rg_to_int[read.get_tag('RG')]

    valid_positions = (original_quals >= minscore)
    qcov = original_quals[valid_positions]
    cycle = bamread_cycle_covariates(read)[valid_positions]
    dinuccov = bamread_dinuc_covariates(read, dinuc_to_int, complement)[valid_positions]

    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] + qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuccov] + positiondeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

"""
given sequences from fastq, get the recalibration vectors
given: names, rawquals, rcorrected, seqs, rgs, seqlen
the vectors are: meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total
"""
def rcorrected_to_vectors(names, rawquals, rcorrected, seqs, rgs, seqlen):
    firstread = np.isin(names,'/1')
    secondread = np.logical_not(firstread)
    pass

def main():
    np.seterr(all = 'raise')
    print(tstamp(), "Starting . . .", file=sys.stderr)
    uncorrfile = "nospace.reads.fq"
    corrfile = "nospace.lighter.fq"
    bamfilename = "only_confident.sorted.recal.bam"
    fastafilename = "../chr1.renamed.fa"
    names, rawquals, rcorrected, seqs, rgs, seqlen = find_rcorrected_sites(uncorrfile, corrfile)

    #rcorrected is true when the bases match the original, false otherwise
    #hence rcorrected == false is where the errors are

    bad_positions = load_positions("variable_sites.txt")
    cachefile = 'cached_recal_errs.npz'
    tablefile = 'only_confident.sorted.recal.txt'
    if os.path.exists(cachefile):
        print(tstamp(), "Loading cached errors . . .", file=sys.stderr)
        loaded = np.load(cachefile)
        gatkcalibratedquals = loaded['gatkcalibratedquals']
        erroneous = loaded['erroneous']
        skips = loaded['skips']
    else:
        gatkcalibratedquals, erroneous, skips = find_errors(bamfilename, fastafilename, bad_positions, names, seqlen)
        np.savez_compressed(cachefile, gatkcalibratedquals = gatkcalibratedquals, erroneous = erroneous, skips = skips)

    #important arrays: names, rawquals, rcorrected, calibquals, gatkcalibratedquals, erroneous, hmmquals

    dinucleotide, dinuc_order = get_dinucleotide(seqs, rawquals)
    unique_rgs = np.unique(rgs)

    rg_to_int = dict(zip(unique_rgs, range(len(unique_rgs))))
    rgs = np.array([rg_to_int[r] for r in rgs], dtype = np.int_)
    rgs = np.broadcast_to(rgs[:,np.newaxis], rgs.shape + (seqlen,)).copy()

    bamfile = pysam.AlignmentFile(bamfilename,"r")
    id_to_pu = {rg['ID'] : rg['PU'] for rg in bamfile.header.as_dict()['RG']}
    unique_pus = [id_to_pu[rg] for rg in unique_rgs]

    # np.set_printoptions(threshold = np.inf)
    # bam_test("only_confident.sorted.recal.bam", tablefile, rg_to_int, unique_pus, dinuc_order, seqlen)
    # quit()

    reversecycle = np.zeros(len(names), dtype = np.bool)
    reversecycle[np.array(list(names.values()), dtype = np.int)] = np.char.endswith(np.array(list(names.keys()), dtype = np.unicode),'/2')

    dq_calibrated = delta_q_recalibrate(rawquals, rgs, dinucleotide, np.logical_not(rcorrected), reversecycle)
    custom_gatk_calibrated = delta_q_recalibrate(rawquals, rgs, dinucleotide, erroneous, reversecycle)
    from_table = table_recalibrate(rawquals, tablefile, unique_pus, dinuc_order, seqlen, reversecycle, rgs, dinucleotide)
    try:
        assert np.array_equal(from_table, gatkcalibratedquals)
    except AssertionError:
        print("from_table shape:", from_table.shape)
        print("gatkcalibratedquals shape:", gatkcalibratedquals.shape)
        print("gatkcalibratedquals[0,:]", gatkcalibratedquals[0,:])
        print("from_table[0,:]", from_table[0,:])
        print("rawquals[0,:]", rawquals[0,:])
        print("foundinplp[0,:]", foundinplp[0,:])
        ne = (from_table != gatkcalibratedquals)
        print("ne.shape",ne.shape)
        print("sum(ne)", np.sum(ne))
        print("sum(from_table.mask != gatkcalibratedquals.mask)", np.sum(from_table.mask != gatkcalibratedquals.mask))
        raise

    #nonsnp = (np.sum(erroneous, axis = 1) > 1) #reads with more than 1 "error" (ie an indel)
    #skips[nonsnp,:] = True

    print(tstamp(), "Skipping", np.sum(skips), "of", skips.size, "(", np.sum(skips)/skips.size ,"%)", "Sites . . .", file=sys.stderr)
    raw = rawquals[~skips]
    gatk = gatkcalibratedquals[~skips]
    dq = dq_calibrated[~skips]
    custom = custom_gatk_calibrated[~skips]
    table = from_table[~skips]
    truth = erroneous[~skips]

    plot_calibration([raw, gatk, dq, custom],
        truth = truth,
        labels = ["Uncalibrated Scores", "GATK Calibration", "KBBQ", "GATK Python Implementation"],
        plotname = 'qualscores.png',
        plottitle = "Substitution Error Calibration")


if __name__ == '__main__':
    main()

