"""
A mish-mash of functions for recalibration. Should probably be renamed in the near future.
"""

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
import kbbq.recaltable as recaltable
import kbbq.plot
import pandas as pd
import datetime
import kbbq.benchmark

def tstamp():
    return '[ ' + datetime.datetime.today().isoformat(' ', 'seconds') + ' ]'

def load_positions(posfile):
    """
    Use with a BED file
    """
    d = dict()
    with open(posfile, 'r') as infh:
        for line in infh:
            # bed format: pos is 0 based and end is 1 based
            chrom, pos, end = line.rstrip().split()
            for i in range(int(pos), int(end)):
                d.setdefault(chrom, list()).append(i)
    return d

def get_var_sites(vcf):
    """
    Use with a VCF file
    """
    vcf = pysam.VariantFile(vcf)
    d = dict()
    for record in vcf:
        for i in range(record.start, record.stop, 1):
            #record.start is 0-based inclusive
            #record.stop is 0-based exclusive
            d.setdefault(record.chrom, list()).append(i)
    return d

def find_corrected_sites(uncorrfile, corrfile):
    print(tstamp(), "Finding corrected sites . . .", file=sys.stderr)
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
    corrected = np.zeros([len(uncorr_reads),seqlen], dtype = np.bool)
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
        corrected[i] = (uncorr_s == corr_s)
    return names, rawquals.copy(), corrected.copy(), seqs.copy(), rgs.copy(), seqlen

def train_regression(rawquals, corrected, tol = 1e-4):
    print(tstamp(), "Doing Logit Regression", file=sys.stderr)
    lr = LR(tol = tol)
    lr = lr.fit(rawquals.flatten().reshape(-1,1), corrected.flatten())
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
            #gatk counts all insertions as aligning to the ref base to the right (maybe?)
            #i think for now we skip if it's skipped on both sides
            skips[readidx:readidx+l] = np.logical_and(subset_variable[refidx-1], subset_variable[refidx])
            readidx = readidx + l
        elif op == 2 or op == 3:
            #deletion in read or N op
            # N is for introns in mRNA
            skips[readidx - 1] = np.logical_or(skips[readidx-1],np.any(subset_variable[refidx: refidx + l]))
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
            raise ValueError("Unrecognized Cigar Operation " + str(op) + " In Read\n" + str(read))
    return readerrors, skips

class RescaledNormal:
    """
    A class to cache the rescaled normal prior used in the bayesian recalibration
    model. 

    Attributes

        * :attr:`maxscore` - max score supported
        * :attr:`prior_dist` - numpy array of prior probability 

    Methods

        * :meth:`prior` - get the prior probability for a given quality score difference

    Most of these attributes are nonsense; the "proper" way to interact
    with the class is via the :attr:`prior_dist` array. If that makes you
    uncomfortable, use the provided accessor function :meth:`prior`.

    Under no circumstances should you attempt to replace any of these attributes,
    it will most likely not have the desired effect. The caching mechanism here 
    only works because the class attributes are immediately instatiated when the class
    is created, so by the time you replace them it won't matter. 
    Reimplement it yourself if you want a different prior.
    """

    oldset = np.seterr(all = 'raise')

    maxscore = 42
    """
    The maximum quality score supported by this class.
    """

    possible_diffs = np.arange(maxscore+1, dtype = np.int_)
    prior_dist = np.zeros(possible_diffs.shape[0], dtype = np.longdouble)
    for i in range(possible_diffs.shape[0]):
        try:
            prior_dist[i] = np.log(.9 * np.exp(-((possible_diffs[i]/.5)**2)/2))
        except FloatingPointError:
            prior_dist[i] = np.NINF
    np.seterr(**oldset)

    def prior(self, difference):
        """
        Return the prior probability for a given difference in quality score.

        :param int difference: The difference in quality score
        :returns: the prior probability
        :rtype: np.longdouble
        """
        return prior_dist[difference]

class Dinucleotide:
    """
    A class to cache dinucleotides and maintain a consistent dinuc -> int
    map throughout the module.
    """

    nucleotides = ['A','T','G','C']
    """
    List of valid nucleotides
    """

    complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    """
    Dictionary for complementing nucleotides
    """

    #dinucs = [i + j for i in nucleotides for j in nucleotides]
    #the above won't work because list comprehensions ignore
    #class scope except the outermost variable. we can use
    # a temporary function, even though it looks bad :(
    dinucs = lambda d:[i + j for i in d for j in d]
    dinucs = dinucs(nucleotides)
    """
    List of valid dinucleotides
    """

    dinuc_to_int = dict(zip(dinucs, range(len(dinucs))))
    """
    Dictionary mapping dinuc -> int
    """

    vectorized_get = np.vectorize(dinuc_to_int.get, otypes = [np.int])

    @classmethod
    def vecget(cls, *args, **kwargs):
        return cls.vectorized_get(*args, **kwargs)

def gatk_delta_q(prior_q, numerrs, numtotal, maxscore = 42):
    """
    Calculate the shift in quality scores from the prior given
    data.

    This is achieved by finding the difference between the maximum
    a posteriori and the prior point estimate of Q.
    """
    assert prior_q.shape == numerrs.shape == numtotal.shape
    possible_q = np.arange(maxscore+1, dtype = np.int)
    diff = np.absolute(np.subtract.outer(possible_q, prior_q).astype(np.int))
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
    assert posterior_q.shape == prior_q.shape
    return posterior_q - prior_q



def p_to_q(p, maxscore = 42):
    q = np.zeros(p.shape, dtype = np.int)
    q[p != 0] = (-10.0*np.log10(p[p != 0])).astype(np.int) #avoid divide by 0
    q[p == 0] = maxscore
    q = np.clip(q, 0, maxscore)
    return q.copy()

def q_to_p(q):
    p = np.array(np.power(10.0,-(q / 10.0)), dtype = np.longdouble, copy = True)
    return p

## Generic covariate functions

def generic_cycle_covariate(sequencelen, secondinpair = False):
    cycle = np.arange(sequencelen)
    if secondinpair:
        cycle = np.negative(cycle + 1)
    return cycle

def generic_dinuc_covariate(sequences, quals, minscore = 6):
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
    dinuccov[...,1:][...,np.logical_not(invalid)] = Dinucleotide.vecget(dinuc[...,np.logical_not(invalid)])
    return dinuccov

## Recalibrating FASTQ reads

def fastq_cycle_covariates(read, secondinpair = False):
    return generic_cycle_covariate(len(read.sequence), secondinpair)

def fastq_dinuc_covariates(read, minscore = 6):
    quals = np.array(read.get_quality_array(), dtype = np.int)
    return generic_dinuc_covariate(np.array(list(read.sequence)), quals, minscore)

def fastq_infer_secondinpair(read):
    namestr = read.name.split(sep='_')[0]
    return namestr[-2:] == '/2'

def fastq_infer_rg(read):
    """
    Infer the read group from appended read information, such as produced by the
    samtools fastq command.

    Requires the read to be formatted such the rg tag is added to the end of the
    read name delimited by a ``_``. Returns the read group tag.
    """
    rgstr = read.name.split(sep='_')[1]
    assert rgstr[0:2] == 'RG'
    return rgstr.split(':')[-1]

def recalibrate_fastq(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg, dinuc_to_int, secondinpair = False, minscore = 6, maxscore = 42):
    qcov = np.array(read.get_quality_array(), dtype = np.int)
    recalibrated_quals = np.array(qcov, copy = True, dtype = np.int)
    valid_positions = (qcov >= minscore)
    qcov = qcov[valid_positions]
    cycle = fastq_cycle_covariates(read, secondinpair)[valid_positions]
    dinuccov = fastq_dinuc_covariates(read, minscore)[valid_positions]
    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] + qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuccov] + positiondeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

## Recalibrate reads from a BAM

def bamread_get_oq(read):
    #TODO: add an assert or logic to ensure OQ is present
    oq = np.array(list(read.get_tag('OQ')), dtype = np.unicode)
    quals = np.array(oq.view(np.uint32) - 33, dtype = np.uint32)
    return quals

def get_rg_to_pu(bamfileobj):
    rg_to_pu = {rg['ID'] : rg['PU'] for rg in bamfileobj.header.as_dict()['RG']}
    return rg_to_pu
