"""
Utilities for recalibrating reads.
"""

from kbbq import compare_reads as utils
from kbbq import recaltable
from kbbq.gatk import applybqsr
import kbbq.read
from kbbq import covariate
import pysam
import numpy as np
import pathlib
import subprocess
import sys
import contextlib

def find_corrected_sites(uncorr_read, corr_read):
    """
    Given a read and a corrected read, return an array of corrected sites.

    :param uncorr_read: Uncorrected read
    :type uncorr_read: :class:`kbbq.read.ReadData`
    :param corr_read: Corrected version of the same read
    :type corr_read: :class:`kbbq.read.ReadData`
    :return: Array of changed sites
    :rtype: :class:`numpy.ndarray` of bool
    """
    assert corr_read.name.startswith(uncorr_read.name)
    return (uncorr_read.seq != corr_read.seq)

def recalibrate_read(read, dqs, minscore = 6):
    """
    Return new qualities given :class:`ReadData` and DQ arrays.
    """
    meanq, globaldeltaq, qscoredeltaq, dinucdeltaq, cycledeltaq = dqs
    rg = read.get_rg_int()
    qcov = read.qual
    recalibrated_quals = np.array(qcov, copy = True)
    valid_positions = (qcov >= minscore)
    qcov = qcov[valid_positions]
    cycle = read.get_cycle_array()[valid_positions]
    dinuc = read.get_dinuc_array()[valid_positions]
    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] +\
        qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuc] +\
        cycledeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

def recalibrate_fastq(fastq, dqs, out, infer_rg = False):
    """
    Recalibrate reads in a FASTQ file given a fastq file name and
    :class:`CovariateData` dqs. Out should be a file-like object.
    """
    # dqs = applybqsr.get_delta_qs_from_covariates(covariates)
    rg = fastq if infer_rg is False else None
    with pysam.FastxFile(uncorr) as fin:
        for fqread in fin:
            r = kbbq.read.ReadData.from_fastq(fqread, rg = rg)
            recalibrated_quals = recalibrate_read(r, dqs)
            fqread.quality = list(recalibrated_quals)
            out.write(str(fqread))

def recalibrate_bam(bam, dqs, out):
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
                kbbq.read.load_rgs_from_bamfile(fin)
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
                kbbq.read.ReadData.load_rgs_from_bamfile(fin)
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
        convert = kbbq.read.ReadData.from_bamread
    elif isinstance(iterable pysam.FastxFile):
        convert = kbbq.read.ReadData.from_fastq
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

