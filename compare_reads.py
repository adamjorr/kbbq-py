#!/usr/bin/env python3
import pysam
import numpy as np
import sklearn
from sklearn.linear_model import LogisticRegression as LR
from sklearn.isotonic import IsotonicRegression as IR
import importlib.util
spec = importlib.util.spec_from_file_location("ek", "/home/ajorr1/bin/jelly/error_kmers.py")
ek = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ek)
import os
import os.path
import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import khmer
import pystan
import scipy.stats
spec = importlib.util.spec_from_file_location("recaltable", "/home/ajorr1/bin/kbbq/recaltable.py")
recaltable = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recaltable)
import pandas as pd

def load_positions(posfile):
    d = dict()
    with open(posfile, 'r') as infh:
        for line in infh:
            chrom, pos = line.rstrip().split()
            d.setdefault(chrom, list()).append(int(pos))
    return d

def find_rcorrected_sites(uncorrfile, corrfile):
    print(ek.tstamp(), "Finding Rcorrected sites . . .", file=sys.stderr)
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
    rawquals = np.zeros([len(uncorr_reads), seqlen], dtype = np.int64)
    rcorrected = np.zeros([len(uncorr_reads),seqlen], dtype = np.bool_)
    seqs = np.zeros(len(uncorr_reads), dtype = 'U' + str(seqlen))
    rgs = np.zeros(len(uncorr_reads), dtype = 'U7')
    for i in range(len(uncorr_reads)):
        names[uncorr_reads[i].name.split(sep='_')[0]] = i
        #print(uncorr_reads[i].name)
        rgs[i] = uncorr_reads[i].name.split(sep='_')[1].split(':')[-1]
        rawquals[i,:] = uncorr_reads[i].get_quality_array()
        seqs[i] = uncorr_reads[i].sequence
        uncorr_s = np.array(list(uncorr_reads[i].sequence), dtype = np.unicode_)
        corr_s = np.array(list(corr_reads[i].sequence), dtype = np.unicode_)
        rcorrected[i] = (uncorr_s == corr_s)
    return names, rawquals.copy(), rcorrected.copy(), seqs.copy(), rgs.copy(), seqlen

def train_regression(rawquals, rcorrected, tol = 1e-4):
    print(ek.tstamp(), "Doing Logit Regression", file=sys.stderr)
    lr = LR(tol = tol)
    lr = lr.fit(rawquals.flatten().reshape(-1,1), rcorrected.flatten())
    return lr

def recalibrate(lr, oldquals):
    print(ek.tstamp(), "Recalibrating Quality Scores . . .", file = sys.stderr)
    shape = oldquals.shape
    q = np.ma.masked_array(oldquals, copy = True)
    notmasked = q[~q.mask]
    newprobs = lr.predict_proba(notmasked.flatten().reshape(-1,1))[:,1]
    newq = np.ma.masked_array(-10.0 * np.ma.log10(newprobs), dtype = np.longdouble)
    newq = np.ma.masked_array(np.rint(newq), dtype=np.int)
    newq = np.clip(newq, 0, 43)
    q.flat[np.logical_not(q.mask.flat)] = newq
    assert q.shape == shape
    return q.copy()

def process_plp(plpfilename, var_pos, names, seqlen, suffix):
    """
    Returns an array of gatk calibrated qualities and actually erroneous sites
    """
    print(ek.tstamp(), "Processing Pileup " + plpfilename + " . . .", file = sys.stderr)
    gatkcalibratedquals = np.zeros([len(names), seqlen], dtype = np.int)
    erroneous = np.zeros([len(names), seqlen], dtype = np.bool_)
    trackingmask = np.zeros([len(names), seqlen], dtype = np.bool_)

    numfinder = re.compile('[\+-](\d+)')
    with open(plpfilename, 'r') as infh:
        varidx = {k : 0 for k in var_pos.keys()}
        for line in infh:
            chrom, pos, refbase, depth, bases, quals, qpos, qname = line.split('\t')
            pos = int(pos)
            if var_pos[chrom][varidx[chrom]] == pos:
                #increment unless it's the last position in the chromosome
                varidx[chrom] = (varidx[chrom] + 1 if varidx[chrom] < len(var_pos[chrom]) - 1 else varidx[chrom])
                continue
            if int(depth) == 0:
                continue
            assert refbase != 'N'
            assert var_pos[chrom][varidx[chrom]] > pos
            bases = re.sub('\^.', '', bases)
            bases = bases.replace('$','')
            bases = bases.replace('<','')
            bases = bases.replace('>','')
            
            match = numfinder.search(bases)
            while match:
                #remove only the specified number of bases inserted or deleted
                bases = bases[:match.start()] + bases[(match.end() + int(match[1])):]
                match = numfinder.search(bases)
            
            assert len(bases) == int(depth)
            assert len(bases) == len(quals)

            bases = np.array(list(bases), dtype = np.unicode_)
            errs = np.array(np.logical_and(bases != '.', bases != ','))

            qpos = np.array(qpos.split(','), dtype = np.int) - 1
            qname = qname.rstrip()
            qname = np.array(qname.split(','), dtype = np.unicode_)
            qname = np.core.defchararray.add(qname, suffix)

            quals = np.array(list(quals), dtype = np.unicode_)
            quals = np.array(quals.view(np.uint32) - 33, dtype = np.uint32)
            ##here
            mask = np.array([names.get(n) for n in qname])
            present = np.flatnonzero(mask)
            mask = mask[present].astype(np.int)
            gatkcalibratedquals[mask,qpos[present]] = quals[present]
            erroneous[mask,qpos[present]] = errs[present]
            trackingmask[mask,qpos[present]] = True

    return gatkcalibratedquals.copy(), erroneous.copy(), trackingmask.copy()

class RescaledNormal:
    oldset = np.seterr(all = 'raise')
    maxscore = 43
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

def v_gatk_delta_q(prior_q, numerrs, numtotal, maxscore = 43):
    assert prior_q.shape == numerrs.shape == numtotal.shape
    possible_q = np.arange(maxscore, dtype = np.int)
    diff = np.absolute(np.subtract.outer(possible_q, prior_q).astype(np.int64))
    #1st dim is possible qs
    prior = RescaledNormal.prior_dist[diff]
    #figure out how to make this work
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

def table_to_vectors(tablefile, rg_order, dinuc_order, seqlen, maxscore = 43):
    #the recal table uses the PU of the read group as the read group entry in the table
    table = recaltable.RecalibrationReport(tablefile)
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

def table_recalibrate(q, tablefile, rg_order, dinuc_order, seqlen, reversecycle, rgs, dinucleotide, maxscore = 43):
    meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total = table_to_vectors(tablefile, rg_order, dinuc_order, seqlen, maxscore)
    globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq = get_delta_qs(meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total)
    qmask = np.ma.getmaskarray(q)
    recal_q = np.zeros(q.shape, dtype = np.longdouble)
    recal_q = np.ma.masked_where(qmask, recal_q)
    recal_q = np.ma.masked_where(q <= 6, recal_q)
    pos = np.broadcast_to(np.arange(q.shape[1]), (q.shape[0], q.shape[1])).copy()
    np.add.at(pos, reversecycle, 1)
    np.negative.at(pos,reversecycle)
    recal_q = meanq[rgs] + globaldeltaq[rgs] + qscoredeltaq[rgs,q] + positiondeltaq[rgs, q, pos] + dinucdeltaq[rgs, q, dinucleotide]
    r_q = np.ma.masked_array(np.rint(np.clip(recal_q,0,maxscore)), dtype = np.int, copy = True)

    return r_q.copy()


def delta_q_recalibrate(q, rgs, dinucleotide, errors, reversecycle, maxscore = 43):
    print(ek.tstamp(), "Getting Covariate Arrays . . .", file=sys.stderr)
    meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total = v_get_covariate_arrays(q, rgs, dinucleotide, errors, reversecycle)
    print(ek.tstamp(), "Finding Delta Q's . . .", file=sys.stderr)
    globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq = get_delta_qs(meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total)
    print(ek.tstamp(), "Recalibrating . . .", file=sys.stderr)
    qmask = np.ma.getmaskarray(q)
    recal_q = np.zeros(q.shape, dtype = np.longdouble)
    #recal_q = np.ma.masked_array(np.zeros(q.shape, dtype = np.longdouble), copy = True) #+ meanq + globaldeltaq
    recal_q = np.ma.masked_where(qmask, recal_q)
    #recal_q[q <= 6] = np.ma.masked
    recal_q = np.ma.masked_where(q <= 6, recal_q)

    #vectorization:
    pos = np.broadcast_to(np.arange(q.shape[1]), (q.shape[0], q.shape[1])).copy()
    np.add.at(pos, reversecycle, 1)
    np.negative.at(pos,reversecycle)
    recal_q = meanq[rgs] + globaldeltaq[rgs] + qscoredeltaq[rgs,q] + positiondeltaq[rgs, q, pos] + dinucdeltaq[rgs, q, dinucleotide]

    #clip and round
    r_q = np.ma.masked_array(np.rint(np.clip(recal_q,0,maxscore)), dtype = np.int, copy = True)

    return r_q.copy()

def get_dinucleotide(seqs, q):
    #[A, T, G, C] -> [A, T, G, C]
    #nucleotides at the beginning of the sequence have an empty string before them
    #we should: ignore any context containing an N, ignore any context at beginning of sequence
    #we also ignore the longest string at the start and end of the sequence that have scores <=2
    print(ek.tstamp(), "Getting dinucleotide context . . .", file=sys.stderr)
    nucs = ['A','T','G','C']
    seqs = seqs.copy() #we may need to alter this
    dinucs = [i + j for i in  nucs for j in nucs]
    dinuc_to_int = dict(zip(dinucs, range(len(dinucs))))
    dinucleotide = generic_dinuc_covariate(seqs.view('U1').reshape((seqs.size, -1)), q, dinuc_to_int)
    return dinucleotide.copy(), dinucs.copy()

def v_get_covariate_arrays(q, rgs, dinucleotide, errors, reversecycle, maxscore = 43, minscore = 6):
    #input arrays are the same dimensions: (numsequences, seqlen)
    #output arrays are dimension (nrgs), (nrgs, q), (nrgs, q, seqlen), or (nrgs, q, 16)
    m = np.ma.getmaskarray(q)
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
    e = np.logical_and(errors, ~m)
    rge = rgs[e]
    qe = q[e]
    valid = np.logical_and(dinucleotide != -1, ~m)
    e_and_valid = np.logical_and(errors, valid)

    np.add.at(expected_errs, rgs[~m], q_to_p(q[~m]))
    np.add.at(rg_errs, rge, 1)
    np.add.at(rg_total, rgs[~m], 1)
    np.add.at(q_errs, (rge, qe), 1)
    np.add.at(q_total, (rgs[~m], q[~m]), 1)
    np.add.at(pos_errs, (rge, qe, pos[e]), 1)
    np.add.at(pos_total, (rgs[~m], q[~m], pos[~m]), 1)
    np.add.at(dinuc_errs, (rgs[e_and_valid], q[e_and_valid], dinucleotide[e_and_valid]), 1)
    np.add.at(dinuc_total, (rgs[valid], q[valid], dinucleotide[valid]), 1)

    meanq = p_to_q(expected_errs / rg_total)

    return meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

def get_covariate_arrays(q, rgs, dinucleotide, errors, maxscore = 43, minscore = 6):
    #these arrays should be the same dimensions: (numsequences, seqlen)
    #total mean Q, quality score, position, dinucleotide context
    #TODO: figure how to do this smartly per read group
    qmask = np.ma.getmaskarray(q.copy())
    qmask = np.logical_or(qmask, q <= minscore)
    qunmasked = np.logical_not(qmask)
    #meanq = p_to_q(np.ma.sum(q_to_p(q)) / np.sum(qunmasked))
    error_and_unmasked = np.logical_and(qunmasked, errors)

    global_errs = np.sum(errors[qunmasked])
    global_total = np.sum(qunmasked)
    q_errs = np.bincount(q[error_and_unmasked], minlength = maxscore+1)
    q_total = np.bincount(q[qunmasked], minlength = maxscore+1)
    #meanq prior is q of expected # errors / total observations
    expected_errs = np.sum(q_to_p(np.arange(maxscore+1)) * q_total)
    meanq = p_to_q(expected_errs / np.sum(q_total))

    #pos_errs should be 2d: pos_errs[qscore, position] = numerrors
    pos_errs = np.zeros([maxscore + 1, 2 * q.shape[1]])
    pos_total = np.zeros([maxscore + 1, 2 * q.shape[1]])
    for i in range(maxscore + 1):
        #pos_errs = np.bincount(q[error_and_unmasked], minlength = maxscore+1) is vectorized equivalent?
        pos_errs[i,] = np.sum(np.logical_and(error_and_unmasked, q == i), axis = 0)
        #pos_total = np.bincount(q[qunmasked], minlength = maxscore + 1) is vectorized equivalent?
        pos_total[i,] = np.sum(np.logical_and(qunmasked, q == i), axis = 0)

    #dinuc_errs is 2d: dinuc_errs[qscore, dinuc] = numerrors
    dinuc_errs = np.zeros([maxscore + 1, 16])
    dinuc_total = np.zeros([maxscore + 1, 16])
    for i in range(maxscore+1):
        e_unmasked_and_q = np.logical_and(error_and_unmasked, q == i)
        unmasked_and_q = np.logical_and(qunmasked, q == i)
        e_valid = np.logical_and(error_and_unmasked, dinucleotide != -1)
        valid = np.logical_and(unmasked_and_q, dinucleotide != -1)
        dinuc_errs[i,] = np.bincount(dinucleotide[e_valid], minlength = 16)
        dinuc_total[i,] = np.bincount(dinucleotide[valid], minlength = 16)

    return meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

def get_delta_qs(meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total, maxscore = 43):
    # shapes are:
    #   [rg]
    #   [rg, q]
    #   [rg, q, covariate]
    #seqlen = pos_total.shape[2] / 2
    nrgs = meanq.shape[0]
    rgdeltaq = v_gatk_delta_q(meanq, rg_errs, rg_total)
    # the qscoredeltaq is 1d, with the index being the quality score
    #qscoredeltaq = np.array([gatk_delta_q(meanq + globaldeltaq, q_errs[i], q_total[i]) for i in range(maxscore + 1)])
    prior1 = np.broadcast_to((meanq + rgdeltaq)[:,np.newaxis], q_total.shape).copy()
    qscoredeltaq = v_gatk_delta_q( prior1 , q_errs, q_total)
    ## positiondeltaq is 2d, first dimension is quality score and second is position
    #positiondeltaq = np.zeros([maxscore+1, seqlen])
    ## dinucdeltaq is 2d, first dimension is quality score and second is nucleotide context
    #dinucdeltaq = np.zeros([maxscore+1, 17])
    prior2 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], pos_total.shape).copy()
    positiondeltaq = v_gatk_delta_q(prior2, pos_errs, pos_total)
    prior3 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], dinuc_total.shape).copy()
    dinucdeltaq = v_gatk_delta_q(prior3, dinuc_errs, dinuc_total)

    #need to add another value of dinuc, for invalid dinuc
    pad = np.zeros((len(dinucdeltaq.shape),2), dtype = np.int_)
    pad[-1,1] = 1 #add a 0 to the last axis
    dinucdq = np.pad(dinucdeltaq, pad_width = pad, mode = 'constant', constant_values = 0)

    return rgdeltaq.copy(), qscoredeltaq.copy(), positiondeltaq.copy(), dinucdq.copy()

def plot_calibration(data, truth, labels, plotname, plottitle = None):
    print(ek.tstamp(), "Making Quality Score Plot . . .", file = sys.stderr)
    assert np.ndim(truth) == 1
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    qualplot, (ax1, ax2) = plt.subplots(2, sharex = True, gridspec_kw = { 'height_ratios' : [2, 1] }, figsize = (7,9))
    #ax1.set_aspect('equal')
    ax2.set_ylim(top = 3e7)
    qualplot.suptitle(plottitle)
    maxscore = 43
    ax1.plot(np.arange(maxscore), 'k:', label = "Perfect")
    try:
        assert np.all([np.array_equal(data[0].shape, data[i].shape) for i in range(len(data))])
    except AssertionError:
        for i in range(len(data)):
            print(labels[i], 'Shape:', data[i].shape)
    maskstack = np.stack([np.ma.getmaskarray(data[i]) for i in range(len(data))] + [np.ma.getmaskarray(truth)])
    allmasks = np.any(maskstack, axis = 0)
    for i in range(len(data)):
        label = labels[i]
        print(ek.tstamp(), "Plotting %s . . ." % (label), file = sys.stderr)
        assert np.ndim(data[i]) == 1
        estimate = np.ma.masked_where(allmasks, data[i], copy = True)
        est_p = q_to_p(estimate)
        unmask = np.logical_not(allmasks)
        try:
            bscore = sklearn.metrics.brier_score_loss(truth[unmask].reshape(-1), est_p[unmask].reshape(-1))
        except ValueError:
            print(est_p)
            print(unmask)
            print(est_p[unmask])
            raise
        numtotal = np.bincount(estimate[unmask].reshape(-1), minlength = (maxscore+1))
        numerrs = np.bincount(estimate[np.logical_and(unmask, truth)].reshape(-1), minlength = len(numtotal)) #not masked and error
        numtotal = np.ma.masked_equal(numtotal, 0)
        p = np.ma.divide(numerrs,numtotal)
        q = p_to_q(p)
        ax1.plot(np.arange(len(q))[np.logical_not(q.mask)], q[np.logical_not(q.mask)], 'o-', alpha = .6, label ="%s, %1.5f" % (label, bscore))
        ax2.plot(np.arange(len(q))[np.logical_not(q.mask)], numtotal[np.logical_not(q.mask)], 'o-', alpha = .6, label = "%s, %1.5f" % (label, bscore))
    plt.xlabel("Predicted Quality Score")
    ax1.set_ylabel("Actual Quality Score")
    ax1.legend(loc = "upper left")
    ax2.set_ylabel("Sample Size")
    qualplot.savefig(plotname)

def p_to_q(p, maxscore = 43):
    q = -10.0*np.ma.log10(p)
    q = np.ma.masked_array(np.rint(q), dtype=np.int, copy = True)
    q = np.clip(q, 0, maxscore)
    return q.copy()

def q_to_p(q):
    p = np.ma.masked_array(np.ma.power(10.0,-(q / 10.0)), dtype = np.longdouble, copy = True)
    return p

def load_pileups(plpfile1, plpfile2, bad_positions, names, seqlen):
    gatkcalibratedquals1, erroneous1, trackingmask1 = process_plp(plpfile1, bad_positions, names, seqlen, "/1")
    gatkcalibratedquals2, erroneous2, trackingmask2 = process_plp(plpfile2, bad_positions, names, seqlen, "/2")
    foundinplp = np.logical_or(trackingmask1,trackingmask2)
    reversecycle = np.ma.masked_where(~foundinplp, trackingmask2)
    gatkcalibratedquals = np.ma.masked_where(~foundinplp, gatkcalibratedquals1 + gatkcalibratedquals2)
    gatkcalibratedquals[gatkcalibratedquals == 2] = np.ma.masked
    erroneous = np.ma.masked_where(~foundinplp, np.logical_or(erroneous1, erroneous2))
    return foundinplp, gatkcalibratedquals, erroneous, reversecycle



def bam_test(bamfile, tablefile, rg_to_int, rg_order, dinuc_order, seqlen, minscore = 6, maxscore = 43):
    print(ek.tstamp(), "Beginning BAM test . . .", file = sys.stderr)
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
    print(ek.tstamp(), "BAM test completed successfully.", file = sys.stderr)



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
    dinuccov[...,1:][...,np.logical_not(invalid)] = vecget(dinuccov[...,1:][...,np.logical_not(invalid)])
    return dinuccov

## Recalibrating FASTQ reads

def fastq_cycle_covariates(read, secondinpair = False):
    return generic_cycle_covariate(len(read.sequence), secondinpair)

def fastq_dinuc_covariates(read, dinuc_to_int, minscore = 6):
    quals = np.array(read.get_quality_array(), dtype = np.int)
    return generic_dinuc_covariate(read.sequence, quals, dinuc_to_int, minscore)

def recalibrate_fastq(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg, dinuc_to_int, secondinpair = False, minscore = 6, maxscore = 43):
    qcov = np.array(read.get_quality_array(), dtype = np.int)
    recalibrated_quals = np.array(qcov, copy = True, dtype = np.int)
    valid_positions = (qcov >= minscore)
    cycle = fastq_cycle_covariates(read)[valid_positions]
    dinuccov = fastq_dinuc_covariates(read, dinuc_to_int, minscore)[valid_positions]
    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] + qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuccov] + positiondeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

## Recalibrate reads from a BAM

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

def recalibrate_bamread(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg_to_int, dinuc_to_int, minscore = 6, maxscore = 43):
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
    print(ek.tstamp(), "Starting . . .", file=sys.stderr)
    uncorrfile = "nospace.reads.fq"
    corrfile = "nospace.reads.cor.fq"
    names, rawquals, rcorrected, seqs, rgs, seqlen = find_rcorrected_sites(uncorrfile, corrfile)
    rawquals = np.ma.masked_equal(rawquals,2) #2 is not a quality score in this data

    #rcorrected is true when the bases match the original, false otherwise
    #hence rcorrected == false is where the errors are

    bad_positions = load_positions("variable_sites.txt")
    plpfile1 = "only_confident.recal.1.plp"
    plpfile2 = "only_confident.recal.2.plp"
    cachefile = 'cached_recal_errs.npz'
    tablefile = 'only_confident.sorted.recal.txt'
    if os.path.exists(cachefile):
        print(ek.tstamp(), "Loading cached errors . . .", file=sys.stderr)
        loaded = np.load(cachefile)
        foundinplp = loaded['foundinplp']
        assert not np.all(~foundinplp)
        gatkcalibratedquals = loaded['gatkcalibratedquals']
        gatkcalibratedquals = np.ma.masked_where(~foundinplp, gatkcalibratedquals)
        gatkcalibratedquals[gatkcalibratedquals == 2] = np.ma.masked
        erroneous = loaded['erroneous']
        erroneous = np.ma.masked_where(~foundinplp, erroneous)
        reversecycle = loaded['reversecycle']
        reversecycle = np.ma.masked_where(~foundinplp, reversecycle)
    else:
        foundinplp, gatkcalibratedquals, erroneous, reversecycle = load_pileups(plpfile1, plpfile2, bad_positions, names, seqlen)
        assert not np.all(~foundinplp)
        np.savez_compressed(cachefile, foundinplp = foundinplp, gatkcalibratedquals = gatkcalibratedquals, erroneous = erroneous, reversecycle = reversecycle)

    #important arrays: names, rawquals, rcorrected, calibquals, gatkcalibratedquals, erroneous, hmmquals

    dinucleotide, dinuc_order = get_dinucleotide(seqs, rawquals)
    unique_rgs = np.unique(rgs)

    rg_to_int = dict(zip(unique_rgs, range(len(unique_rgs))))
    rgs = np.array([rg_to_int[r] for r in rgs], dtype = np.int_)
    rgs = np.broadcast_to(rgs[:,np.newaxis], rgs.shape + (seqlen,)).copy()

    bamfile = pysam.AlignmentFile("only_confident.sorted.recal.bam","r")
    id_to_pu = {rg['ID'] : rg['PU'] for rg in bamfile.header.as_dict()['RG']}
    unique_pus = [id_to_pu[rg] for rg in unique_rgs]

    # np.set_printoptions(threshold = np.inf)
    # bam_test("only_confident.sorted.recal.bam", tablefile, rg_to_int, unique_pus, dinuc_order, seqlen)
    # quit()

    dq_calibrated = delta_q_recalibrate(rawquals.copy(), rgs, dinucleotide, np.logical_not(rcorrected), reversecycle)
    custom_gatk_calibrated = delta_q_recalibrate(rawquals.copy(), rgs, dinucleotide, erroneous, reversecycle)
    from_table = table_recalibrate(rawquals.copy(), tablefile, unique_pus, dinuc_order, seqlen, reversecycle, rgs, dinucleotide)

    plot_calibration([rawquals.flatten(), gatkcalibratedquals.flatten(), dq_calibrated.flatten(), custom_gatk_calibrated.flatten(), from_table.flatten()],
        truth = erroneous.flatten(),
        labels = ["Uncalibrated Scores", "GATK Calibration", "KBBQ", "GATK Python Implementation", "From Table"],
        plotname = 'qualscores.png',
        plottitle = "Comparison of Calibration Methods")


if __name__ == '__main__':
    main()

