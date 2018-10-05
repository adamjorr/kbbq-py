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
    for i in range(len(uncorr_reads)):
        names[uncorr_reads[i].name] = i
        rawquals[i,:] = uncorr_reads[i].get_quality_array()
        seqs[i] = uncorr_reads[i].sequence
        uncorr_s = np.array(list(uncorr_reads[i].sequence), dtype = np.unicode_)
        corr_s = np.array(list(corr_reads[i].sequence), dtype = np.unicode_)
        rcorrected[i] = (uncorr_s == corr_s)
    return names, rawquals.copy(), rcorrected.copy(), seqs.copy(), seqlen

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

def plot_calibration(data, truth, labels, plotname, plottitle = None):
    print(ek.tstamp(), "Making Quality Score Plot . . .", file = sys.stderr)
    assert np.ndim(truth) == 1
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    qualplot = plt.figure()
    ax = qualplot.add_subplot(111)
    qualplot.suptitle(plottitle)
    maxscore = 43
    plt.plot(np.arange(maxscore), 'k:', label = "Perfect")
    allmasks = np.stack([data[i].mask for i in range(len(data))].append(truth.mask))
    allmasks = np.any(allmasks, axis = 0)
    for i in range(len(data)):
        label = labels[i]
        print(ek.tstamp(), "Plotting %s . . ." % (label), file = sys.stderr)
        assert np.ndim(data[i]) == 1
        estimate = np.ma.masked_array(data[i])
        estimate[allmasks] = np.ma.masked
        est_p = np.ma.masked_array(np.ma.power(10.0,-(estimate / 10.0)), dtype = np.longdouble)
        unmask = np.logical_not(est_p.mask)
        bscore = sklearn.metrics.brier_score_loss(truth[unmask].reshape(-1), est_p[unmask].reshape(-1))
        numtotal = np.bincount(estimate[unmask].reshape(-1), minlength = (maxscore+1))
        numerrs = np.bincount(estimate[unmask][truth[unmask]], minlength = len(numtotal)) #not masked and error
        numtotal = np.ma.masked_equal(numtotal, 0)
        p = np.ma.divide(numerrs,numtotal)
        q = -10.0*np.ma.log10(p)
        q = np.ma.masked_array(np.rint(q), dtype=np.int)
        q = np.clip(q, 0, maxscore)
        plt.plot(np.arange(len(q))[np.logical_not(q.mask)], q[np.logical_not(q.mask)], 'o-', label ="%s, %1.3f" % (label, bscore))
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Actual Quality Score")
    plt.legend(loc = "upper left")
    qualplot.savefig(plotname)

def p_to_q(p, maxscore = 43):
    q = -10.0*np.ma.log10(p)
    q = np.ma.masked_array(np.rint(q), dtype=np.int)
    q = np.clip(q, 0, maxscore)
    return q

def q_to_p(q):
    p = np.ma.masked_array(np.ma.power(10.0,-(q / 10.0)), dtype = np.longdouble)
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

def pystan_model(modelfile, L, seqlen, ksize, alpha, beta, a, raw_p):
    os.environ['STAN_NUM_THREADS'] = "32"
    extra_compile_args = ['-pthreads', '-DSTAN_THREADS']
    sm = pystan.StanModel(file = modelfile, model_name = 'kbbq', extra_compile_args = extra_compile_args)
    datadict = {'L' : L, 'seqlen' : seqlen, 'ksize' : ksize, 'alpha': alpha, 'beta': beta, 'a' : a}
    estimates = sm.optimizing(data = datadict,
        sample_file = 'samples.csv',
        init = lambda chain_id = None: {'e' : raw_p, 'lambda' : [30, 1300]},
        verbose = True)
    return estimates['e'], estimates['kerr'], estimates['lambda']

def run_stan_model(modelfile, seqlen, ksize, rawquals, abundances):
    raw_p = q_to_p(rawquals.copy())
    adjustable = np.logical_not(np.any(raw_p.mask, axis = 1))
    alpha, beta, _, _ = scipy.stats.beta.fit(raw_p[adjustable,:].data.astype(np.float), floc = 0, fscale = 1) ##needs log1p which isn't implemented for longdouble type.
    a = abundances[adjustable,:]
    L = a.shape[0]
    e, kerr, l = pystan_model(modelfile, L, seqlen, ksize, alpha, beta, a, raw_p[adjustable,:])
    raw_p[adjustable,:] = e
    recalibrated = p_to_q(raw_p.copy())
    recalibrated[~adjustable,:] = np.ma.masked
    return recalibrated

def main():
    np.seterr(all = 'raise')
    print(ek.tstamp(), "Starting . . .", file=sys.stderr)
    uncorrfile = "reads.fq"
    corrfile = "nospace.reads.cor.fq"
    khmerfile = "counts.khmer"
    names, rawquals, rcorrected, seqs, seqlen = find_rcorrected_sites(uncorrfile, corrfile)
    rawquals = np.ma.masked_equal(rawquals,2) #2 is not a quality score in this data

    lr = train_regression(rawquals[~rawquals.mask], np.logical_not(rcorrected[~rawquals.mask]))
    #rcorrected is true when the bases match the original, false otherwise
    #hence rcorrected == false is where the errors are
    calibquals = recalibrate(lr, rawquals)

    #qtostr = np.arange(43, dtype = np.uint32) + 33
    #qtostr = qtostr.view('U1')

    bad_positions = load_positions("variable_sites.txt")
    plpfile1 = "only_confident.1.plp"
    plpfile2 = "only_confident.2.plp"
    gatkcalibratedquals1, erroneous1, trackingmask1 = process_plp(plpfile1, bad_positions, names, seqlen, "/1")
    gatkcalibratedquals2, erroneous2, trackingmask2 = process_plp(plpfile2, bad_positions, names, seqlen, "/2")
    foundinplp = np.logical_or(trackingmask1,trackingmask2)
    gatkcalibratedquals = np.ma.masked_where(~foundinplp, gatkcalibratedquals1 + gatkcalibratedquals2)
    gatkcalibratedquals[gatkcalibratedquals == 2] = np.ma.masked
    erroneous = np.ma.masked_where(~foundinplp, np.logical_or(erroneous1, erroneous2))

    ## This is a cheesy function to run the test without a lot of boilerplate
    #readnames, readquals = ek.do_recalibration()
    #hmmquals = np.ma.masked_array(np.zeros([len(names), seqlen], dtype = np.int64))
    #for i in range(len(readquals)):
    #    idx = names.get(readnames[i])
    #    if idx is None:
    #        continue
    #    else:
    #        hmmquals[idx,:] = readquals[i]
    #notinhmm = np.all(hmmquals == 0, axis = 1)
    #hmmquals[notinhmm,:] = np.ma.masked

    #important arrays: names, rawquals, rcorrected, calibquals, gatkcalibratedquals, erroneous, hmmquals

    ##do a regression on the hmm
    #tunmask = np.logical_not(np.logical_or(erroneous.mask, hmmquals.mask))
    #lr2 = train_regression(hmmquals[~hmmquals.mask], np.logical_not(rcorrected[~hmmquals.mask]))
    #lr2recal = recalibrate(lr2, hmmquals)

    ##do an isotonic regression
    #ir = IR(out_of_bounds = 'clip')
    #rawp = np.ma.power(10.0, -(rawquals / 10.0))
    #isop = np.ma.masked_array(rawquals, dtype = np.longdouble)
    #isop[~isop.mask] = ir.fit_transform(rawp[~rawp.mask].data, np.logical_not(rcorrected[~rawquals.mask].data))
    #isoquals = -10.0*np.ma.log10(isop)
    #isoquals = np.ma.masked_array(np.rint(isoquals), dtype=np.int)
    #isoquals = np.clip(isoquals, 0, 43)

    ##try to roll your own naive bayes classifier
    ##naivebayes = khmer_nb(seqs, rawquals, erroneous, seqlen, khmerfile)

    abundances, ksize, nkmers = get_abundances(seqs, seqlen, khmerfile, savefile = 'abundances.txt')
    ##def pystan_model(modelfile, L, seqlen, ksize, prior_e, a)
    ## def run_stan_model(modelfile, seqlen, ksize, rawquals, abundances):
    stan_calibrated = run_stan_model('/home/ajorr1/bin/kbbq/kbbq.stan', seqlen, ksize, rawquals, abundances)

    plot_calibration([rawquals.flatten(), gatkcalibratedquals.flatten(), calibquals.flatten(), rcorrected.flatten()*rawquals.flatten(), stan_calibrated.flatten()],
        truth = erroneous.flatten(),
        labels = ["Uncalibrated Scores", "GATK Calibration", "KBBQ - Logit Regression", "Rcorrector", "KBBQ - MAP"],
        plotname = 'qualscores.png',
        plottitle = "Comparison of Calibration Methods")


if __name__ == '__main__':
    main()

