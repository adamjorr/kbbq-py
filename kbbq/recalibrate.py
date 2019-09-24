"""
Utilities for recalibrating reads.
"""

import kbbq
from kbbq import compare_reads as utils
from kbbq import recaltable
import kbbq.gatk.applybqsr
import kbbq.read
import kbbq.covariate
import kbbq.bloom
import khmer
import pysam
import numpy as np
import pathlib
import subprocess
import sys
import contextlib
import itertools

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
    meanq, globaldeltaq, qscoredeltaq, cycledeltaq, dinucdeltaq = dqs
    rg = read.get_rg_int()
    qcov = read.qual
    recalibrated_quals = np.array(qcov, copy = True, dtype = np.int)
    valid_positions = (qcov >= minscore)
    qcov = qcov[valid_positions]
    cycle = read.get_cycle_array()[valid_positions]
    dinuc = read.get_dinucleotide_array()[valid_positions]
    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] +\
        qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuc] +\
        cycledeltaq[rg, qcov, cycle]).astype(np.int)
    recalibrated_quals = np.clip(recalibrated_quals, 0, utils.RescaledNormal.maxscore)
    return recalibrated_quals

def recalibrate_fastq(fastq, dqs, out, infer_rg = False):
    """
    Recalibrate reads in a FASTQ file given a fastq file name and
    :class:`CovariateData` dqs. Out should be a file-like object.
    """
    rg = fastq if infer_rg is False else None
    with pysam.FastxFile(uncorr) as fin:
        for fqread in fin:
            r = kbbq.read.ReadData.from_fastq(fqread, rg = rg)
            recalibrated_quals = recalibrate_read(r, dqs)
            fqread.quality = list(recalibrated_quals)
            out.write(str(fqread))

def find_covariates(read_sources):
    """
    Consume reads from the given list of iterables and load
    them into a :class:`kbbq.covariate.CovariateData` object.
    """
    data = kbbq.covariate.CovariateData()
    for r in itertools.chain.from_iterables(read_sources):
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
                kbbq.read.ReadData.load_rgs_from_bamfile(fin)
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
def open_outputs(files, output, bams):
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
                header = fin.header.to_dict()
                pgids = [x.get('ID','') for x in header['PG']]
                suffixes = [int(x.split('.',1)[-1]) for x in pgids if x.startswith("kbbq")]
                headerid = 'kbbq' if suffixes == [] else 'kbbq.' + str(max(suffixes) + 1)
                header['PG'] = header['PG'] + [{'ID' : headerid, 'PN' : 'kbbq', 'CL' : ' '.join(sys.argv), 'VN' : kbbq.__version__}]
            fout = pysam.AlignmentFile(str(o), mode = 'wb', header = header)
        else:
            fout = open(str(o), mode = 'w')
        opened_outputs.append(fout)
    try:
        yield opened_outputs
    finally:
        for o in opened_outputs:
            o.close()

def yield_reads(iterable, *args, **kwargs):
    """
    Return a generator of (ReadData, original_read) pairs.
    """
    if isinstance(iterable, pysam.AlignmentFile):
        convert = kbbq.read.ReadData.from_bamread
    elif isinstance(iterable, pysam.FastxFile):
        convert = kbbq.read.ReadData.from_fastq
    else:
        raise ValueError("Unknown iterable type {}".format(type(iterable)))
    for i in iterable:
        yield convert(i, *args, **kwargs), i

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
            generators.append(yield_reads(f, use_oq = use_oq))
        else: #fastq file
            if infer_rg == False:
                rg = str(i)
            else:
                rg = None
            generators.append(yield_reads(f, rg = rg))
    try:
        yield generators
    finally:
        for f in opened_files:
            f.close()

def recalibrate(files, output, infer_rg = False, use_oq = False, set_oq = False, ksize = 32, memory = '2G', alpha = .1, gatkreport = None):
    #make these options later
    if output == []:
        output = [str(i.with_name(i.stem + '.kbbq' + i.suffix)) for i in [pathlib.Path(f) for f in files ]]
    if len(files) != len(output):
        raise ValueError('One output should be specified for each input.')
    validate_files(files)
    bams = load_headers_from_bams(files)

    if gatkreport is None or not pathlib.Path(gatkreport).is_file():
        #gatkreport not provided or the provided report doesn't exist
        utils.print_info("Loading hash")
        graph = kbbq.bloom.create_empty_nodegraph(ksize = ksize, max_mem = memory)
        with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData generators
            #1st pass: load hash
            for read, original in itertools.chain.from_iterable(allreaddata): #a single ReadData generator
                kbbq.bloom.count_read(read, graph, sampling_rate = alpha)

        fpr = khmer.calc_expected_collisions(graph, force = False, max_false_pos = .15)
        utils.print_info("False positive rate:", str(fpr))
        p_added = kbbq.bloom.p_kmer_added(sampling_rate = alpha, graph = graph)
        utils.print_info("Probability any k-mer was added:", str(p_added))
        thresholds = kbbq.bloom.calculate_thresholds(p_added, graph.ksize())
        utils.print_info("Error thresholds:", thresholds)
        covariates = kbbq.covariate.CovariateData()
        utils.print_info("Finding trusted k-mers")
        trustgraph = kbbq.bloom.create_empty_nodegraph(ksize = ksize, max_mem = memory)
        with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData generators
            #pass 1.5: find trusted kmers
            for read, original in itertools.chain.from_iterable(allreaddata): #a single ReadData generator
                errors = kbbq.bloom.infer_read_errors(read, graph, thresholds)
                #don't trust bad quality
                errors[read.qual <= 6] = True
                read.errors = errors
                #find k-size blocks of non-errors and add to the trusted graph
                kbbq.bloom.add_trusted_kmers(read, trustgraph)

        utils.print_info("Finding errors and building model")
        with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData generators
            #2nd pass: find errors + build model
            for read, original in itertools.chain.from_iterable(allreaddata): #a single ReadData generator
                read.errors = kbbq.bloom.infer_errors_from_trusted_kmers(read, trustgraph)
                covariates.consume_read(read)

        dqs = kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)
        #TODO: if gatkreport doesn't exist, save the model to it
    else:
        #gatkreport provided and exists
        #TODO: update to use new DQ / Covariate API
        #TODO: this won't work with a fastq because we don't know the RGs yet!!
        #we need to get the RGs from the table and load them rather than the other way
        #around
        utils.print_info("Loading model from", str(gatkreport))
        meanq, *recalvecs = kbbq.gatk.applybqsr.table_to_vectors(kbbq.recaltable.RecalibrationReport.fromfile(gatkreport), rg_order = list(kbbq.read.ReadData.rg_to_pu.values()))
        dqs = kbbq.gatk.applybqsr.ModelDQs(meanq, *kbbq.gatk.applybqsr.get_delta_qs(meanq, *recalvecs))

    utils.print_info("Recalibrating reads")
    with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata, \
        open_outputs(files, output, bams) as opened_outputs: #a list of ReadData generators
        #3rd pass: recalibrate
        for i, o in zip(allreaddata, opened_outputs):
            for read, original in i:
                recalibrated_quals = recalibrate_read(read, dqs)
                if isinstance(original, pysam.AlignedSegment):
                    if set_oq:
                        original.set_tag('OQ',
                            pysam.array_to_qualitystring(original.query_qualities))
                    if original.is_reverse:
                        #the read class considers the read as fwd, so the recalibrated
                        #quals will be flipped if the read is reversed.
                        recalibrated_quals = np.flip(recalibrated_quals)
                    original.query_qualities = list(recalibrated_quals)
                    o.write(original)
                elif isinstance(original, pysam.FastxRecord):
                    original.quality = utils.nparray_to_qualitystring(recalibrated_quals)
                    o.write(str(original) + '\n')
                else:
                    raise ValueError("Unknown read type {}".format(type(original)))
    utils.print_info("Done")
