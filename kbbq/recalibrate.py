"""
Utilities for recalibrating reads.
"""

from kbbq import compare_reads as utils
from kbbq import recaltable
from kbbq.gatk import applybqsr
import pysam
import numpy as np
import pathlib
import subprocess
import sys
import contextlib

#TODO: probably add 2 classes: a ReadData class and a CovariateData class

def find_corrected_sites(uncorr_read, corr_read):
    """
    Given a fastq read and a corrected fastq read, return an array of corrected sites.
    """
    assert corr_read.name.startswith(uncorr_read.name)
    uncorr_seq = np.array(list(uncorr_read.sequence), dtype = np.unicode)
    corr_seq = np.array(list(corr_read.sequence), dtype = np.unicode)
    return (uncorr_seq != corr_seq)

def fastq_to_covariate_arrays(fastq, infer_rg = False, minscore = 6, maxscore = 42):
    """
    TODO:We really really really need a class to keep track of the covariate arrays
    and manage their size. This is ridiculous. This function is ugly as sin.
    """
    #initialize dict for keeping track of read groups
    if infer_rg is False:
        rgfun = lambda x: 0 #since we currently only support 1 fastq at a time
    else:
        rgfun = utils.fastq_infer_rg
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
            pos = utils.fastq_cycle_covariates(uncorr_read, utils.fastq_infer_secondinpair(uncorr_read))
            dinucleotide = utils.fastq_dinuc_covariates(uncorr_read, minscore)

            skips = np.zeros(seqlen, dtype = np.bool)
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
            np.add.at(expected_errs, rgv, utils.q_to_p(qv))
            np.add.at(rg_errs, rge, 1)
            np.add.at(rg_total, rgv, 1)
            np.add.at(q_errs, (rge, qe), 1)
            np.add.at(q_total, (rgv, qv), 1)
            np.add.at(pos_errs, (rge, qe, pos[e_and_valid]), 1)
            np.add.at(pos_total, (rgv, qv, pos[valid]), 1)
            np.add.at(dinuc_errs, (rgs[e_and_dvalid], q[e_and_dvalid], dinucleotide[e_and_dvalid]), 1)
            np.add.at(dinuc_total, (rgs[dinuc_valid], q[dinuc_valid], dinucleotide[dinuc_valid]), 1)
    meanq = utils.p_to_q(expected_errs / rg_total)
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
        rgfun = utils.fastq_infer_rg
    rg_to_int = dict()
    nrgs = 0

    meanq, *vectors = fastq_to_covariate_arrays(fastq, infer_rg)
    dqs = applybqsr.get_delta_qs(meanq, *vectors)
    with pysam.FastxFile(fastq[0]) as fin:
        for read in fin:
            rg = rgfun(read)
            rgint = rg_to_int.get(rg)
            if rgint is None:
                rgint = nrgs
                rg_to_int[rg] = rgint
                nrgs = nrgs + 1
            recalibrated_quals = utils.recalibrate_fastq(read, meanq, *dqs,
                rg = rgint, dinuc_to_int = utils.Dinucleotide.dinuc_to_int,
                secondinpair = utils.fastq_infer_secondinpair(read))
            strquals = ''.join((recalibrated_quals + 33).astype(np.uint32).view('U1'))
            print('@' + read.name)
            print(read.sequence)
            print('+')
            print(strquals)

def recalibrate_bam(read, output, original):
    """
    Recalibrate ReadData read and output to filehandle output.
    """


def find_covariates(read_sources):
    """
    Consume reads from the given list of iterables and load
    them into a :class:`kbbq.covariate.CovariateData` object.
    """
    data = covariate.CovariateData()
    for i in read_sources:
        for r in i:
            data.consume_read(r)
    return data

def opens_as_bam(path):
    """
    Attempt to open the path as a BAM file. Return True if no error is raised.
    """
    try:
        with pysam.AlignmentFile(str(path)) as fin:
            pass
        return True
    except ValueError:
        return False

def load_headers_from_bams(inputs):
    """
    Load RG information from the headers of each input bam. Return which inputs are bams.
    """
    bams = []
    for i in inputs:
        if opens_as_bam(i):
            bams.append(True)
            with pysam.AlignmentFile(str(i)) as fin:
                read.load_rgs_from_bamfile(fin)
        else:
            bams.append(False)
    return bams

def validate_files(files):
    """
    Ensures input files are regular files.
    """
    for i in files:
        path = pathlib.Path(i)
            if not path.is_file():
                raise ValueError('Given path {} does not exist or is not a regular file'.format(path))

@contextlib.contextmanager
def open_outputs(files, output, bams, infer_rg, use_oq, set_oq):
    """
    Initialize output files ensuring BAM headers are handled.

    This function acts as a context manager. The outputs
    will automatically be closed once they go out of context.

    See https://docs.python.org/3/library/contextlib.html for
    more information.
    """
    opened_outputs = []
    for i, o, b in zip(files, output, bams):
        if b:
            with pysam.AlignmentFile(i) as fin:
                header = fin.header
                headerid = max([x.get('ID',0) for x in header['PG']]) + 1
                header['PG'] = header['PG'] + [{'ID' : headerid, 'PN' : 'kbbq', 'CL' : ' '.join(sys.argv)}]
            fout = pysam.AlignmentFile(str(o), mode = 'wb', header = header)
        else:
            fout = pysam.FastxFile(str(o), mode = 'w')
        opened_outputs.append(fout)
    try:
        yield opened_outputs
    finally:
        for o in opened_outputs:
            o.close()

def yield_reads(iterable, *args, **kwargs):
    """
    Return a generator of ReadData objects.
    """
    if isinstance(iterable, pysam.AlignmentFile):
        convert = read.ReadData.from_bamread
    elif isinstance(iterable pysam.FastxFile):
        convert = read.ReadData.from_fastq
    else:
        raise ValueError("Unknown iterable type {}".format(type(iterable)))
    for i in iterable:
        yield convert(i, *args, **kwargs)

@contextlib.contextmanager
def generate_reads_from_files(files, bams, infer_rg = False, use_oq = False):
    """
    Return a list of ReadData generators from opening each file.

    This function acts as a context manager. The files
    will automatically be closed once they go out of context.

    See https://docs.python.org/3/library/contextlib.html for
    more information.
    """
    opened_files = [pysam.AlignmentFile(i) if b else pysam.FastxFile(i) for i,b in zip(files,bams)]
    generators = []
    for f,b,i in zip(opened_files,bams,files):
        if b: #bam file
            generators.append(yield_reads(f, use_oq))
        else: #fastq file
            if infer_rg = False:
                rg = i
            else:
                rg = None
            generators.append(yield_reads(f, rg = rg))
    # for r in fin:
    #     if set_oq is True:
    #         r.set_tag(tag = 'OQ', value = ''.join([chr(q + 33) for q in r.query_qualities]), value_type = 'Z')
    #     yield read.ReadData.from_bamread(read, use_oq), fout, r
    try:
        yield generators
    finally:
        for f in opened_files:
            f.close()

def recalibrate(files, output, cmd, infer_rg = False, use_oq = False, set_oq = False, gatkreport = None):
    if gatkreport is not None:
        raise NotImplementedError('GATKreport reading / creation is not yet supported.')
    if output == []:
        output = [ str(i.with_name(i.stem + '.kbbq' + i.suffix)) for pathlib.Path(i) in files ]
    if len(files) != len(output):
        raise ValueError('One output should be specified for each input.')
    validate_files(files)
    bams = load_headers_from_bams(files)
    with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData iterators
        #1st pass: load hash

        #2nd pass: find errors + build model

        #3rd pass: recalibrate


    for r, output, original in generate_reads_from_files(files, output, cmd, infer_rg, use_oq, set_oq):
        #ReadData, output, and original read (AlignedSegment or FastxProxy)
        # 3 passes:
        # 1st pass: load hash
        # 2nd pass: find errors + build model
        # 3rd pass: recalibrate
    with open_outputs(files, output, bams, infer_rg, use_oq, set_oq) as opened_outputs:
        #do things with open outputs here

