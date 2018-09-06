#!/usr/bin/env python3
import pysam
import numpy as np
from sklearn.linear_model import LogisticRegression as LR
from sklearn.isotonic import IsotonicRegression as IR
import importlib.util
spec = importlib.util.spec_from_file_location("ek", "/home/ajorr1/bin/jelly/error_kmers.py")
ek = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ek)
import os.path
import sys
import re

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

    names = np.zeros(len(uncorr_reads), dtype = np.unicode_)
    seqlen = len(uncorr_reads[0].get_quality_array())
    rawquals = np.zeros([len(uncorr_reads), seqlen], dtype = np.int64)
    rcorrected = np.zeros([len(uncorr_reads),seqlen], dtype = np.bool_)
    for i in range(len(uncorr_reads)):
        names[i] = uncorr_reads[i].name
        rawquals[i,:] = uncorr_reads[i].get_quality_array()
        uncorr_s = np.array(list(uncorr_reads[i].sequence), dtype = np.unicode_)
        corr_s = np.array(list(corr_reads[i].sequence), dtype = np.unicode_)
        rcorrected[i] = (uncorr_s == corr_s)
    return names, rawquals, rcorrected, seqlen

def train_regression(rawquals, rcorrected, tol = 1e-4):
    print(ek.tstamp(), "Doing Logit Regression", file=sys.stderr)
    lr = LR(tol = tol)
    lr.fit(rawquals.flatten().reshape(-1,1), rcorrected.flatten())
    return lr

def recalibrate(lr, quals):
    print(ek.tstamp(), "Recalibrating Quality Scores . . .", file = sys.stderr)
    shape = quals.shape
    probs = np.array(lr.predict_proba(quals.flatten().reshape(-1,1))[:,1].reshape(shape), dtype = np.longdouble)
    q = np.array(-10.0 * np.log10(probs), dtype = np.longdouble)
    quals = np.array(np.rint(q), dtype=np.int)
    quals = np.clip(quals, 0, 43)
    return quals

def process_plp(plpfilename, var_pos, names, seqlen, suffix):
    """
    Returns an array of gatk calibrated qualities and actually erroneous sites
    """
    print(ek.tstamp(), "Processing Pileup " + plpfilename + " . . .", file = sys.stderr)
    gatkcalibratedquals = np.zeros([len(names), seqlen], dtype = np.int)
    erroneous = np.zeros([len(names), seqlen], dtype = np.bool_)

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
            mask = np.isin(names, qname)
            gatkcalibratedquals[mask,qpos] = quals
            erroneous[mask,qpos] = errs

    return gatkcalibratedquals.copy(), erroneous.copy()

def plot_qual_scores(numerrors, numtotal, plotname, plottitle = None):
    print(ek.tstamp(), "Making Base Quality Score Plot . . .", file=sys.stderr)
    plottitle = (plottitle if plottitle is not None else plotname)
    numtotal = np.ma.masked_array(numtotal, mask = (numtotal == 0))
    #numtotal[numtotal == 0] = np.nan #don't divide by 0
    p = numerrors/numtotal
    #p[p == 0] = 1e-4 #1e-4 is the largest quality score, 40
    q = -10.0*np.ma.log10(p)
    q = np.ma.masked_array(np.rint(q), dtype=np.int)
    q = np.clip(q, 0, 42)
    x = np.arange(len(p))
    mse = np.mean(np.square(x - q))
    
    sns.set()
    qualplot = plt.figure()
    ax = qualplot.add_subplot(111)
    qualplot.suptitle(plottitle)
    plt.plot(x)
    plt.plot(np.arange(len(q))[np.logical_not(q.mask)], q[np.logical_not(q.mask)] , 'o-')
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Actual Quality Score")
    plt.legend(labels = ["Perfect","Estimated"], loc = "upper left")
    plt.text(0.5,0.01, "Mean squared error: " + str(mse),
        horizontalalignment = 'center', verticalalignment = 'bottom',
        transform = ax.transAxes)
    qualplot.savefig(plotname)

def main():
    np.seterr(all = 'raise')
    print(ek.tstamp(), "Starting . . .", file=sys.stderr)
    uncorrfile = "reads.fq"
    corrfile = "nospace.reads.cor.fq"
    names, rawquals, rcorrected, seqlen = find_rcorrected_sites(uncorrfile, corrfile)
    lr = train_regression(rawquals, rcorrected)
    calibquals = recalibrate(lr, rawquals)

    #qtostr = np.arange(43, dtype = np.uint32) + 33
    #qtostr = qtostr.view('U1')

    bad_positions = load_positions("variable_sites.txt")
    plpfile1 = "only_confident.1.plp"
    plpfile2 = "only_confident.2.plp"
    gatkcalibratedquals1, erroneous1 = process_plp(plpfile1, bad_positions, names, seqlen, "/1")
    gatkcalibratedquals2, erroneous2 = process_plp(plpfile2, bad_positions, names, seqlen, "/2")
    gatkcalibratedquals = (gatkcalibratedquals1.flatten() + gatkcalibratedquals2.flatten()).reshape(gatkcalibratedquals1.shape)
    erroneous = np.logical_or(erroneous1.flatten(), erroneous2.flatten()).reshape(erroneous1.shape)

    #important arrays: names, rawquals, rcorrected, calibquals, gatkcalibratedquals, erroneous
    print(ek.tstamp(), "Tallying . . .", file=sys.stderr)
    numerrors = np.zeros(42, dtype = np.uint64)
    numtotal = np.zeros(42, dtype = np.uint64)
    np.add.at(numerrors, rawquals.flatten()[erroneous.flatten()], 1)
    np.add.at(numtotal, rawquals.flatten(), 1)

    caliberrs = np.zeros(42, dtype = np.uint64)
    calibtotal = np.zeros(42, dtype = np.uint64)
    np.add.at(caliberrs, calibquals.flatten()[erroneous.flatten()], 1)
    np.add.at(calibtotal, calibquals.flatten(), 1)

    ek.plot_qual_scores(numerrs, numtotal, "qualscores.png", "Raw Reads")
    ek.plot_qual_scores(caliberrs, calibtotal, "calibrated.png", "After Calibration")
    assert np.all(numerrs <= numtotal)
    assert np.all(caliberrs <= calibtotal)
    assert np.all(numtotal == calibtotal)

if __name__ == '__main__':
    main()

