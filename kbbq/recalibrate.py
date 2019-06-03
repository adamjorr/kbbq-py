"""
Utilities for recalibrating reads.
"""

from kbbq import compare_reads
from kbbq import recaltable
import pysam
import numpy as np

#TODO: probably add 2 classes: a ReadData class and a CovariateData class

def find_corrected_sites(uncorr_read, corr_read):
    """
    Given a fastq read and a corrected fastq read, return an array of corrected sites.
    """
    assert corr_read.name.startswith(uncorr_read.name)
    uncorr_seq = np.array(list(uncorr_read.sequence), dtype = np.unicode)
    corr_seq = np.array(list(corr_read.sequence), dtype = np.unicode)
    return (uncorr_seq != corr_seq)

def get_fq_skips(read):
    """
    This function exists just to override it in testing
    """
    seqlen = len(list(read.sequence))
    return np.zeros(seqlen, dtype = np.bool)

def fastq_to_covariate_arrays(fastq, infer_rg = False, minscore = 6, maxscore = 42):
    """
    TODO:We really really really need a class to keep track of the covariate arrays
    and manage their size. This is ridiculous. This function is ugly as sin.
    """
    #initialize dict for keeping track of read groups
    if infer_rg is False:
        rgfun = lambda x: 0 #since we currently only support 1 fastq at a time
    else:
        rgfun = compare_reads.fastq_infer_rg
    rg_to_int = dict()
    nrgs = 0
    seqlen = 0

    qshape = [nrgs, maxscore + 1]
    posshape = [nrgs, maxscore + 1, 2 * seqlen]
    dinucshape = [nrgs, maxscore + 1, 16]

    #initialize covariate arrays
    #these will be minimally initialized and will need to be
    #expanded; however, there should be at least 1 rg (dim 0)
    #and 2 position covariates (if each read is 1 bp long, dim 2).
    meanq = np.zeros(nrgs, dtype = np.int)
    expected_errs = np.zeros(nrgs, dtype = np.longdouble)

    rg_errs = np.zeros(nrgs, dtype = np.int)
    rg_total = np.zeros(nrgs, dtype = np.int)
    q_errs = np.zeros(qshape, dtype = np.int)
    q_total = np.zeros(qshape, dtype = np.int)
    pos_errs = np.zeros(posshape, dtype = np.int)
    pos_total = np.zeros(posshape, dtype = np.int)
    dinuc_errs = np.zeros(dinucshape, dtype = np.int)
    dinuc_total = np.zeros(dinucshape, dtype = np.int)

    with pysam.FastxFile(fastq[0]) as uncorr_in, pysam.FastxFile(fastq[1]) as corr_in:
        for uncorr_read, corr_read in zip(uncorr_in, corr_in):
            #get rg (as an int) and update rg_to_int dict
            rg = rgfun(uncorr_read)
            rgint = rg_to_int.get(rg)
            if rgint is None:
                rgint = nrgs
                rg_to_int[rg] = rgint
                nrgs = nrgs + 1
                #resize all the arrays since we have a new rg
                qshape[0] = qshape[0] + 1
                posshape[0] = posshape[0] + 1
                dinucshape[0] = dinucshape[0] + 1

                meanq = np.append(meanq, 0)
                expected_errs = np.append(expected_errs, 0)
                rg_errs = np.append(rg_errs, 0)
                rg_total = np.append(rg_total, 0)

                q_errs = np.append(q_errs, [np.zeros(qshape[1:], dtype = np.int)], axis = 0)
                q_total = np.append(q_total, [np.zeros(qshape[1:], dtype = np.int)], axis = 0)
                pos_errs = np.append(pos_errs, [np.zeros(posshape[1:], dtype = np.int)], axis = 0)
                pos_total = np.append(pos_total, [np.zeros(posshape[1:], dtype = np.int)], axis = 0)
                dinuc_errs = np.append(dinuc_errs, [np.zeros(dinucshape[1:], dtype = np.int)], axis = 0)
                dinuc_total = np.append(dinuc_total, [np.zeros(dinucshape[1:], dtype = np.int)], axis = 0)
            readlen = len(list(uncorr_read.sequence))
            if readlen > seqlen:
                seqlen = readlen
                padding = 2 * seqlen - posshape[2]
                posshape[2] = 2 * seqlen
                pos_errs = np.append(pos_errs, np.zeros(posshape[0:-1] + [padding], dtype = np.int), axis = 2)
                pos_total = np.append(pos_total, np.zeros(posshape[0:-1] + [padding], dtype = np.int), axis = 2)
            # get covariate values
            rgs = np.zeros(seqlen, dtype = np.int)
            rgs[:] = rgint
            errors = find_corrected_sites(uncorr_read, corr_read)
            q = np.array(uncorr_read.get_quality_array(), dtype = np.int)
            pos = compare_reads.fastq_cycle_covariates(uncorr_read, compare_reads.fastq_infer_secondinpair(uncorr_read))
            dinucleotide = compare_reads.fastq_dinuc_covariates(uncorr_read, compare_reads.Dinucleotide.dinuc_to_int, minscore)

            skips = get_fq_skips(uncorr_read) #this fn exists specifically to be overriden for testing
            skips[q < minscore] = True
            valid = ~skips
            dinuc_valid = np.logical_and(dinucleotide != -1, valid)
            e_and_valid = np.logical_and(errors, valid)
            e_and_dvalid = np.logical_and(errors, dinuc_valid)

            #print(uncorr_read.name, q) # errors.shape, valid.shape, e_and_valid.shape, np.sum(errors), np.sum(valid), np.sum(e_and_valid))

            rge = rgs[e_and_valid]
            rgv = rgs[valid]
            qe = q[e_and_valid]
            qv = q[valid]

            #tally it all up
            np.add.at(expected_errs, rgv, compare_reads.q_to_p(qv))
            np.add.at(rg_errs, rge, 1)
            np.add.at(rg_total, rgv, 1)
            np.add.at(q_errs, (rge, qe), 1)
            np.add.at(q_total, (rgv, qv), 1)
            np.add.at(pos_errs, (rge, qe, pos[e_and_valid]), 1)
            np.add.at(pos_total, (rgv, qv, pos[valid]), 1)
            np.add.at(dinuc_errs, (rgs[e_and_dvalid], q[e_and_dvalid], dinucleotide[e_and_dvalid]), 1)
            np.add.at(dinuc_total, (rgs[dinuc_valid], q[dinuc_valid], dinucleotide[dinuc_valid]), 1)
    meanq = compare_reads.p_to_q(expected_errs / rg_total)
    return meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

def recalibrate_fastq(fastq, infer_rg = False):
    """
    Recalibrate a FASTQ file given a list containing 1) a fastq and 2) a corrected fastq

    TODO: get a better solution for getting all the readgroups out of a fastq
    I think the best thing is going to be using a class variable for a ReadData
    class or something similar. This is going to be absolutely necessary for us
    to eventually support intermixed fastq and bam inputs
    """
    if infer_rg is False:
        rgfun = lambda x: 0 #since we currently only support 1 fastq at a time
    else:
        rgfun = compare_reads.fastq_infer_rg
    rg_to_int = dict()
    nrgs = 0

    meanq, *vectors = fastq_to_covariate_arrays(fastq, infer_rg)
    dqs = compare_reads.get_delta_qs(meanq, *vectors)
    with pysam.FastxFile(fastq[0]) as fin:
        for read in fin:
            rg = rgfun(read)
            rgint = rg_to_int.get(rg)
            if rgint is None:
                rgint = nrgs
                rg_to_int[rg] = rgint
                nrgs = nrgs + 1
            recalibrated_quals = compare_reads.recalibrate_fastq(read, meanq, *dqs,
                rg = rgint, dinuc_to_int = compare_reads.Dinucleotide.dinuc_to_int,
                secondinpair = compare_reads.fastq_infer_secondinpair(read))
            strquals = ''.join((recalibrated_quals + 33).astype(np.uint32).view('U1'))
            print('@' + read.name)
            print(read.sequence)
            print('+')
            print(strquals)

def recalibrate_bam(bam, use_oq = False, set_oq = False):
    """
    Not yet implemented.
    """
    raise NotImplementedError('Recalibrating a bam is not yet implemented. \
        Try converting your BAM to a FASTQ file with the samtools fastq command. \
        We welcome pull requests if you\'d like to work on this feature.')

def recalibrate(bam, fastq, infer_rg = False, use_oq = False, set_oq = False, gatkreport = None):
    if gatkreport is not None:
        raise NotImplementedError('GATKreport reading / creation is not yet supported.')
    elif bam is not None:
        recalibrate_bam(bam, use_oq, set_oq)
    elif fastq is not None:
        recalibrate_fastq(fastq, infer_rg = infer_rg)
