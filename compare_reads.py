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
    diff = np.absolute(np.subtract.outer(possible_q, prior_q))
    #1st dim is possible qs
    prior = RescaledNormal.prior_dist[diff]
    #figure out how to make this work
    #maybe broadcast possible_q at the beginning and use regular subtract ?
    broadcast_errs = np.broadcast_to(numerrs, possible_q.shape + numerrs.shape)
    broadcast_tot = np.broadcast_to(numtotal, possible_q.shape + numtotal.shape)
    p = q_to_p(possible_q).astype(np.float)
    while len(p.shape) < len(broadcast_tot.shape):
        p = np.expand_dims(p, -1)
    broadcast_p = np.broadcast_to( p, broadcast_tot.shape )
    loglike = scipy.stats.binom.logpmf(broadcast_errs+1, broadcast_tot+2, broadcast_p)
    #loglike should now have same dims as prior
    assert loglike.shape == prior.shape
    posterior = prior + loglike
    #this doesn't seem right, check it
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

def gatk_delta_q(prior_q, numerrs, numtotal, maxscore = 43):
    possible_q = np.arange(maxscore, dtype = np.int)
    diff = np.array(np.clip(np.absolute(np.rint(prior_q - possible_q)), 0, maxscore))
    # this is a rescaled normal distribution
    # this underflows if diff >= 20
    prior_dist = np.zeros(diff.shape)
    for i in range(diff.shape[0]):
        try:
            prior_dist[i] = np.log(.9 * np.exp(-((diff[i]/.5)**2)/2))
        except FloatingPointError:
            prior_dist[i] = np.NINF
    #prior_dist = scipy.stats.norm.logpdf(diff, scale = .5)
    #smooth by adding 1 error and 1 nonerror
    #in gatk, this is done by rounding up and adding one
    q_likelihood = scipy.stats.binom.logpmf(numerrs + 1, numtotal+2, q_to_p(possible_q).astype(np.float))
    posterior_q = np.argmax(prior_dist + q_likelihood)
    return posterior_q - prior_q

def delta_q_recalibrate(q, rgs, dinucleotide, errors, maxscore = 43):
    print(ek.tstamp(), "Getting Covariate Arrays . . .", file=sys.stderr)
    meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total = v_get_covariate_arrays(q, rgs, dinucleotide, errors)
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
    posrange = np.arange(recal_q.shape[1])
    recal_q = meanq[rgs] + globaldeltaq[rgs] + qscoredeltaq[rgs,q] + positiondeltaq[rgs, q, posrange] + dinucdeltaq[rgs, q, dinucleotide]

    #clip and round
    r_q = np.ma.masked_array(np.rint(np.clip(recal_q,0,maxscore)), dtype = np.int, copy = True)

    #quantization
    #quantizer = np.arange(maxscore + 1)
    #quantizer[0:3] = 3 #collapse 0, 1, and 2 to 3
    #quantizer[32:] = 32 #collapse everything above 32 to 32
    #recal_q = quantizer[recal_q]

    return r_q.copy()

def get_dinucleotide(q, seqs, seqlen, minq = 2):
    #[A, T, G, C] -> [A, T, G, C]
    #nucleotides at the beginning of the sequence have an empty string before them
    #we should: ignore any context containing an N, ignore any context at beginning of sequence
    #we also ignore the longest string at the start and end of the sequence that have scores <=2
    print(ek.tstamp(), "Getting dinucleotide context . . .", file=sys.stderr)
    nucs = ['A','T','G','C']
    seqs = seqs.copy() #we may need to alter this
    dinucs = [i + j for i in  nucs for j in nucs]
    dinuc_to_int = dict(zip(dinucs, range(len(dinucs))))
    dinucleotide = np.zeros([seqs.shape[0], seqlen], dtype = np.int)
    for i in range(seqs.shape[0]):
        currentseq = seqs[i]
        for j in range(seqlen):
            if q[i,j] <= minq:
                currentseq[j] = 'N'
            else:
                break
        for j in reversed(range(seqlen)):
            if q[i,j] <= minq:
                currentseq[j] = 'N'
            else:
                break

        dinucleotide[i,0] = -1
        for j in range(1, seqlen):
            bases = currentseq[j-1:j+1]
            if not 'N' in bases:
                dinucleotide[i,j] = dinuc_to_int[bases]
            else:
                dinucleotide[i,j] = -1
    return dinucleotide

def v_get_covariate_arrays(q, rgs, dinucleotide, errors, maxscore = 43, minscore = 6):
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
    pos_errs = np.zeros((nrgs, maxscore + 1, seqlen), dtype = np.int_)
    pos_total = np.zeros((nrgs, maxscore + 1, seqlen), dtype = np.int_)
    dinuc_errs = np.zeros((nrgs, maxscore + 1, 16), dtype = np.int_)
    dinuc_total = np.zeros((nrgs, maxscore + 1, 16), dtype = np.int_)
    pos = np.broadcast_to(np.arange(seqlen), (q.shape[0], seqlen))

    #these will be reused a lot; cache them here
    rge = rgs[errors]
    qe = q[errors]
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
    pos_errs = np.zeros([maxscore + 1, q.shape[1]])
    pos_total = np.zeros([maxscore + 1, q.shape[1]])
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
    seqlen = pos_total.shape[2]
    nrgs = meanq.shape[0]
    rgdeltaq = v_gatk_delta_q(meanq, rg_errs, rg_total)
    # the qscoredeltaq is 1d, with the index being the quality score
    #qscoredeltaq = np.array([gatk_delta_q(meanq + globaldeltaq, q_errs[i], q_total[i]) for i in range(maxscore + 1)])
    prior1 = np.broadcast_to((meanq + rgdeltaq)[:,np.newaxis], q_total.shape)
    qscoredeltaq = v_gatk_delta_q( prior1 , q_errs, q_total)
    ## positiondeltaq is 2d, first dimension is quality score and second is position
    #positiondeltaq = np.zeros([maxscore+1, seqlen])
    ## dinucdeltaq is 2d, first dimension is quality score and second is nucleotide context
    #dinucdeltaq = np.zeros([maxscore+1, 17])
    prior2 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], pos_total.shape)
    positiondeltaq = v_gatk_delta_q(prior2, pos_errs, pos_total)
    prior3 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], dinuc_total.shape)
    dinucdeltaq = v_gatk_delta_q(prior3, dinuc_errs, dinuc_total)

    #need to add another value of dinuc, for invalid dinuc
    pad = np.zeros((len(dinucdeltaq.shape),2), dtype = np.int_)
    pad[-1,1] = 1 #add a 0 to the last axis
    dinucdq = np.pad(dinucdeltaq, pad_width = pad, mode = 'constant', constant_values = 0)

    #TODO: could be another layer where position is dependent on qreported while dinuc is dependent on qreported and position
    #TODO: seems like its just a RG problem where my data ACTUALLY has like 30
    #for i in range(maxscore + 1):
    #    for j in range(seqlen):
    #        positiondeltaq[i,j] = gatk_delta_q(meanq + globaldeltaq + qscoredeltaq[i], pos_errs[i,j], pos_total[i,j])
    #    for j in range(16):
    #        dinucdeltaq[i,j] = gatk_delta_q(meanq + globaldeltaq + qscoredeltaq[i], dinuc_errs[i,j], dinuc_total[i,j])
    #    dinucdeltaq[i,-1] = 0 #no contribution for invalid context
    return rgdeltaq.copy(), qscoredeltaq.copy(), positiondeltaq.copy(), dinucdq.copy()

def plot_calibration(data, truth, labels, plotname, plottitle = None):
    print(ek.tstamp(), "Making Quality Score Plot . . .", file = sys.stderr)
    assert np.ndim(truth) == 1
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    qualplot, (ax1, ax2) = plt.subplots(2, sharex = True, gridspec_kw = { 'height_ratios' : [4, 1] }, figsize = (7,9))
    #ax1.set_aspect('equal')
    ax2.set_ylim(top = 1e7)
    qualplot.suptitle(plottitle)
    maxscore = 43
    ax1.plot(np.arange(maxscore), 'k:', label = "Perfect")
    maskstack = np.stack([np.ma.getmaskarray(data[i]) for i in range(len(data))] + [np.ma.getmaskarray(truth)])
    allmasks = np.any(maskstack, axis = 0)
    for i in range(len(data)):
        label = labels[i]
        print(ek.tstamp(), "Plotting %s . . ." % (label), file = sys.stderr)
        assert np.ndim(data[i]) == 1
        estimate = np.ma.masked_array(data[i], copy = True)
        estimate[allmasks] = np.ma.masked
        est_p = np.ma.masked_array(np.ma.power(10.0,-(estimate / 10.0)), dtype = np.longdouble, copy = True)
        unmask = np.logical_not(allmasks)
        try:
            bscore = sklearn.metrics.brier_score_loss(truth[unmask].reshape(-1), est_p[unmask].reshape(-1))
        except ValueError:
            print(est_p)
            print(unmask)
            print(est_p[unmask])
            raise
        numtotal = np.bincount(estimate[unmask].reshape(-1), minlength = (maxscore+1))
        numerrs = np.bincount(estimate[unmask][truth[unmask]], minlength = len(numtotal)) #not masked and error
        numtotal = np.ma.masked_equal(numtotal, 0)
        p = np.ma.divide(numerrs,numtotal)
        q = -10.0*np.ma.log10(p)
        q = np.ma.masked_array(np.rint(q), dtype=np.int)
        q = np.clip(q, 0, maxscore)
        ax1.plot(np.arange(len(q))[np.logical_not(q.mask)], q[np.logical_not(q.mask)], 'o-', label ="%s, %1.5f" % (label, bscore))
        ax2.plot(np.arange(len(q))[np.logical_not(q.mask)], numtotal[np.logical_not(q.mask)], 'o-', label = "%s, %1.5f" % (label, bscore))
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

def get_abundances(seqs, seqlen, khmerfile, savefile = None):
    args = khmer.khmer_args.build_counting_args().parse_args()
    alltable = khmer.khmer_args.create_countgraph(args)
    alltable.load(khmerfile)
    ksize = alltable.ksize()
    nkmers = seqlen - ksize + 1
    if savefile is not None and os.path.exists(savefile):
        print(ek.tstamp(), "Loading K-mer Abundances", file = sys.stderr)
        abundances = np.loadtxt(savefile, dtype = np.int64)
    else:
        print(ek.tstamp(), "Finding K-mer Abundances", file = sys.stderr)
        getter = alltable.get_kmer_counts
        abundances = np.zeros([len(seqs), nkmers], dtype = np.int64)
        for i in range(len(seqs)):
            abundances[i,:] = getter(seqs[i])
        if savefile is not None:
            np.savetxt(savefile, abundances, fmt = '%d')
    return abundances, ksize, nkmers

def khmer_nb(seqs, rawquals, erroneous, seqlen, khmerfile):
    print(ek.tstamp(), "Naively trying a Naive Bayes Classifier...", file = sys.stderr)
    abundances, ksize, nkmers = get_abundances(seqs, seqlen, khmerfile, savefile = 'abundances.txt')
    abundances = np.ma.masked_array(abundances)
    erroneous_mers = np.ma.masked_array(np.zeros([len(seqs), nkmers], dtype = np.bool_))
    for i in range(len(seqs)):
        for j in range(nkmers):
            erroneous_bases = erroneous[i,j:j+ksize]
            erroneous_mers[i,j] = np.any(erroneous_bases)
            if np.ma.getmask(erroneous_mers[i,j]) is not np.ma.nomask and np.any(erroneous_bases.mask):
                erroneous_mers[i,j] = np.ma.masked
                abundances[i,j] = np.ma.masked
    noterr_or_masked = np.logical_not(np.logical_or(erroneous_mers.mask, erroneous_mers.data))
    err_and_notmasked = np.logical_and(erroneous_mers.data, np.logical_not(erroneous_mers.mask))
    erroneous_abundances = np.bincount(abundances[noterr_or_masked].flatten())
    nonerroneous_abundances = np.bincount(abundances[err_and_notmasked].flatten())

    #make the arrays the same length
    if len(erroneous_abundances) > len(nonerroneous_abundances):
        nonerroneous_abundances.resize(erroneous_abundances.shape)
    elif len(nonerroneous_abundances) > len(erroneous_abundances):
        erronous_abundances.resize(nonerroneous_abundances.shape)

    p_a_given_e = np.array(erroneous_abundances / np.sum(erroneous_abundances), dtype = np.longdouble)
    p_a_given_note = np.array(nonerroneous_abundances / np.sum(nonerroneous_abundances), dtype = np.longdouble)
    #P(E | O) = P(O | E) * P(E) / P(O) = sum(P(O | S) * P(S | E)) * P(E) / P(O)
    p_s_given_e = np.array([0,.128,0,.872])
    p_s_given_note = np.array([.424,0,.006,.57])
    #p_s_given_e = s_given_e / np.sum(s_given_e)
    #p_s_given_note = s_given_note / np.sum(s_given_note)
    actual_s_given_e = np.zeros(4, dtype = np.int)
    actual_s_given_note = np.zeros(4, dtype = np.int)
    oldpe = q_to_p(rawquals.copy())
    newpe = oldpe.copy()
    p_e = 0.0108804719256
    print(ek.tstamp(), "Classifying and Recalibrating", file = sys.stderr)
    for i in range(len(seqs)):
        #just the middles for now
        for j in range(nkmers - ksize - 1):
            allidxs = np.array([j, j+1])
            a = abundances[i,allidxs] # dim 2,
            a_given_e = p_a_given_e[a] # dim 2,
            a_given_note = p_a_given_note[a] # dim 2,
            o = np.stack([a_given_note, a_given_e])
            p_o_given_s = np.outer(o[:,0],o[:,1]).flatten() #first state transition
            p_o_given_e = np.sum(p_o_given_s * p_s_given_e)
            p_o_given_note = np.sum(p_o_given_s * p_s_given_note)
            #p_e = oldpe[i,j+ksize]
            p_o = p_o_given_e * p_e + p_o_given_note * (1.0 - p_e)
            p_e_given_o = p_o_given_e * p_e / p_o
            newpe[i,j + ksize] = p_e_given_o

            if np.ma.getmask(oldpe[i,j+ksize]) is not np.ma.nomask and np.ma.getmask(oldpe[i,j+ksize]):
                newpe[i, j+ksize] = np.ma.masked
                actuals = erroneous_mers[i, allidxs]
                actual_nonerrors = np.logical_not(actuals)
                actual_errors = actuals
                actual_o = np.stack([actual_nonerrors, actual_errors])
                actuals_outer = np.outer(actual_o[:,0],actual_o[:,1]).flatten()
                if erroneous[i,j + ksize]:
                    actual_s_given_e = actual_s_given_e + actuals_outer
                else:
                    actual_s_given_note = actual_s_given_note + actuals_outer
    actual_s_given_e = actual_s_given_e / np.sum(actual_s_given_e)
    actual_s_given_note = actual_s_given_note / np.sum(actual_s_given_note)
    # Modeled P(S|E) [0.   0.   0.   0.   0.   0.   0.25 0.25 0.   0.   0.   0.   0.   0. 0.25 0.25]                                            
    # Actual P(S|E) [0.         0.         0.         0.         0.         0. 0.02765633 0.1004291  0.         0.         0.         0. 0.         0.         0.12954866 0.74236591]            
    # Modeled P(S|NOT E) [0.11111111 0.11111111 0.         0.11111111 0.         0. 0.         0.         0.11111111 0.11111111 0.         0.11111111 0.11111111 0.11111111 0.         0.11111111]
    # Actual P(S|NOT E) [3.17953062e-01 2.12157531e-03 0.00000000e+00 1.04366234e-01 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00 4.98386869e-03 2.48924327e-05 0.00000000e+00 1.72464193e-03 2.11434305e-01 1.59378846e-03 0.00000000e+00 3.55797633e-01]
    print("Modeled P(S|E)",p_s_given_e)
    print("Actual P(S|E)",actual_s_given_e)
    print("Modeled P(S|NOT E)", p_s_given_note)
    print("Actual P(S|NOT E)", actual_s_given_note)
    print("Modeled P(E)",p_e)
    print("Actual P(E)",np.sum(erroneous)/np.sum(~erroneous.mask))
    newquals = p_to_q(newpe)
    return newquals

def pystan_model(modelfile, L, seqlen, ksize, alpha, beta, a, raw_p, q):
    os.environ['STAN_NUM_THREADS'] = "32"
    extra_compile_args = ['-pthread', '-DSTAN_THREADS']
    sm = pystan.StanModel(file = modelfile, model_name = 'kbbq', extra_compile_args = extra_compile_args)
    print(ek.tstamp(), "Model Compiled. Starting optimization . . .", file=sys.stderr)
    datadict = {'L' : L, 'seqlen' : seqlen, 'ksize' : ksize, 'a' : a, 'q' : q, 'alpha_init' : alpha, 'beta_init' : beta}
    estimates = sm.optimizing(data = datadict,
        sample_file = 'samples.csv',
        init = lambda chain_id = None: {'e' : raw_p.copy(), 'kerr' : np.zeros((L, seqlen - ksize + 1)), 'lambda' : np.array([30, 1300]), 'alpha' : alpha, 'beta' : beta},
        verbose = True)
    return estimates['e'], estimates['kerr'], estimates['lambda'], estimates['alpha'], estimates['beta']

def run_stan_model(modelfile, seqlen, ksize, rawquals, abundances):
    raw_p = q_to_p(rawquals.copy())
    adjustable = np.logical_not(np.any(raw_p.mask, axis = 1))
    init_alpha, init_beta, _, _ = scipy.stats.beta.fit(raw_p[adjustable,:].data.astype(np.float), floc = 0, fscale = 1) ##needs log1p which isn't implemented for longdouble type.
    a = abundances[adjustable,:]
    L = a.shape[0]
    init_p = raw_p[adjustable,:]
    assert init_p.shape[0] == L
    assert a.shape[0] == L
    assert init_p.shape[1] == seqlen
    assert a.shape[1] == seqlen - ksize + 1
    e, kerr, l, alpha, beta = pystan_model(modelfile, L, seqlen, ksize, init_alpha, init_beta, a, init_p, rawquals[adjustable,:])
    print("Init alpha:", init_alpha)
    print("Init beta:", init_beta)
    print("Fit alpha:", alpha)
    print("Fit beta:", beta)
    print("Init phi:", init_alpha / (init_alpha + init_beta))
    print("Init gamma:", (init_alpha + init_beta))
    print("Sum(perr):", sum(init_p))
    print("Number of elements:", init_p.size)
    print("Mean(perr):", np.mean(init_p))
    print("Var(perr):", np.var(init_p))
    raw_p[adjustable,:] = e
    recalibrated = p_to_q(raw_p.copy())
    recalibrated[~adjustable,:] = np.ma.masked
    return recalibrated

def load_pileups(plpfile1, plpfile2, bad_positions, names, seqlen):
    gatkcalibratedquals1, erroneous1, trackingmask1 = process_plp(plpfile1, bad_positions, names, seqlen, "/1")
    gatkcalibratedquals2, erroneous2, trackingmask2 = process_plp(plpfile2, bad_positions, names, seqlen, "/2")
    foundinplp = np.logical_or(trackingmask1,trackingmask2)
    gatkcalibratedquals = np.ma.masked_where(~foundinplp, gatkcalibratedquals1 + gatkcalibratedquals2)
    gatkcalibratedquals[gatkcalibratedquals == 2] = np.ma.masked
    erroneous = np.ma.masked_where(~foundinplp, np.logical_or(erroneous1, erroneous2))
    return foundinplp, gatkcalibratedquals, erroneous

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
    if os.path.exists(cachefile):
        print(ek.tstamp(), "Loading cached errors . . .", file=sys.stderr)
        loaded = np.load(cachefile)
        foundinplp = loaded['foundinplp']
        assert not np.all(~foundinplp)
        #print("foundinplp.shape:",foundinplp.shape)
        #print("sum(foundinplp):",np.sum(foundinplp))
        gatkcalibratedquals = loaded['gatkcalibratedquals']
        gatkcalibratedquals = np.ma.masked_where(~foundinplp, gatkcalibratedquals)
        gatkcalibratedquals[gatkcalibratedquals == 2] = np.ma.masked
        #print("gatkcalibratedquals.shape:",gatkcalibratedquals.shape)
        #print("sum(np.ma.getmaskarray(gatkcalibratedquals)):",np.sum(np.ma.getmaskarray(gatkcalibratedquals)))
        erroneous = loaded['erroneous']
        erroneous = np.ma.masked_where(~foundinplp, erroneous)
        #print("erroneous.shape:",erroneous.shape)
        #print("sum(np.ma.getmaskarray(erroneous)):",np.sum(np.ma.getmaskarray(gatkcalibratedquals)))
    else:
        foundinplp, gatkcalibratedquals, erroneous = load_pileups(plpfile1, plpfile2, bad_positions, names, seqlen)
        assert not np.all(~foundinplp)
        np.savez_compressed(cachefile, foundinplp = foundinplp, gatkcalibratedquals = gatkcalibratedquals, erroneous = erroneous)

    #important arrays: names, rawquals, rcorrected, calibquals, gatkcalibratedquals, erroneous, hmmquals

    #abundances, ksize, nkmers = get_abundances(seqs, seqlen, khmerfile, savefile = 'abundances.txt')
    ##def pystan_model(modelfile, L, seqlen, ksize, prior_e, a)
    ## def run_stan_model(modelfile, seqlen, ksize, rawquals, abundances):
    #stan_calibrated = run_stan_model('/home/ajorr1/bin/kbbq/kbbq.stan', seqlen, ksize, rawquals, abundances)

    dinucleotide = get_dinucleotide(rawquals, seqs, seqlen)
    unique_rgs = np.unique(rgs)
    #for i in range(unique_rgs.shape[0]):
    #    print(ek.tstamp(), "Processing RG:",unique_rgs[i], "(", i+1, "of", unique_rgs.shape[0], ")",". . .", file=sys.stderr)
    #    thisrg = np.zeros(rawquals.shape, dtype = np.bool_)
    #    thisrg[rgs == unique_rgs[i],] = True
    #    #d_q_recalibrate requires arrays of the right shape and returns them
    #    #so we have to mask the input and only assign the unmasked values
    #    rg_quals = np.ma.masked_where(np.logical_not(thisrg), rawquals)
    #    dq_calibrated[thisrg] = delta_q_recalibrate(rg_quals, rgs, dinucleotide, np.logical_not(rcorrected))[thisrg]
    #    custom_gatk_calibrated[thisrg] = delta_q_recalibrate(rg_quals, rgs, dinucleotide, erroneous)[thisrg]

    rg_to_int = dict(zip(unique_rgs, range(len(unique_rgs))))
    rgs = np.array([rg_to_int[r] for r in rgs], dtype = np.int_)
    rgs = np.broadcast_to(rgs[:,np.newaxis], rgs.shape + (seqlen,))

    dq_calibrated = delta_q_recalibrate(rawquals.copy(), rgs, dinucleotide, np.logical_not(rcorrected))
    custom_gatk_calibrated = delta_q_recalibrate(rawquals.copy(), rgs, dinucleotide, erroneous)

    plot_calibration([rawquals.flatten(), gatkcalibratedquals.flatten(), rcorrected.flatten()*rawquals.flatten(), dq_calibrated.flatten(), custom_gatk_calibrated.flatten()],
        truth = erroneous.flatten(),
        labels = ["Uncalibrated Scores", "GATK Calibration", "Rcorrector", "KBBQ", "GATK Python Implementation"],
        plotname = 'qualscores.png',
        plottitle = "Comparison of Calibration Methods")


if __name__ == '__main__':
    main()

