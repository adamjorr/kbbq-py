"""
Utilities for recalibrating reads.
"""

import kbbq
from kbbq import compare_reads as utils
from kbbq import recaltable
import kbbq.gatk.applybqsr
import kbbq.gatk.bqsr
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
import gzip

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
    recalibrated_quals[valid_positions] = (
        meanq[rg] +\
        globaldeltaq[rg] +\
        qscoredeltaq[rg, qcov] +\
        dinucdeltaq[rg, qcov, dinuc] +\
        cycledeltaq[rg, qcov, cycle]).astype(np.int)
    recalibrated_quals = np.clip(recalibrated_quals, 0, utils.RescaledNormal.maxscore)
    return recalibrated_quals

def recalibrate_fastq(fastq, dqs, out, infer_rg = False):
    """
    Recalibrate reads in a FASTQ file given a fastq file name and
    :class:`CovariateData` dqs. Out should be a file-like object.
    """
    rg = fastq if infer_rg is False else None
    with pysam.FastxFile(fastq) as fin:
        for fqread in fin:
            r = kbbq.read.ReadData.from_fastq(fqread, rg = rg)
            recalibrated_quals = recalibrate_read(r, dqs)
            fqread.quality = utils.nparray_to_qualitystring(recalibrated_quals)
            out.write(str(fqread) + '\n')

def find_covariates(read_sources):
    """
    Consume reads from the given list of iterables and load
    them into a :class:`kbbq.covariate.CovariateData` object.
    """
    data = kbbq.covariate.CovariateData()
    for r in itertools.chain.from_iterable(read_sources):
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

    See :py:mod:`contextlib` for more information.
    """
    opened_outputs = []
    for i, o, b in zip(files, output, bams):
        if b:
            with pysam.AlignmentFile(i) as fin:
                header = fin.header.to_dict()
                oldpg = header.get('PG',[])
                pgids = [x.get('ID','') for x in oldpg]
                kbbqids = [x for x in pgids if x.startswith("kbbq")]
                suffixes = []
                if len(kbbqids) > 0: #multiple IDs
                    #one should have a suffix (but maybe not)
                    for kid in kbbqids:
                        suffix = 0
                        try:
                            suffix = int(kid.split('.',1)[-1])
                        except ValueError:
                            pass
                        suffixes.append(suffix)
                headerid = 'kbbq' if suffixes == [] else 'kbbq.' + str(max(suffixes) + 1)
                header['PG'] = oldpg + [{'ID' : headerid, 'PN' : 'kbbq', 'CL' : ' '.join(sys.argv), 'VN' : kbbq.__version__}]
            fout = pysam.AlignmentFile(str(o), mode = 'wb', header = header)
        else:
            if pathlib.Path(o).suffix == '.gz':
                fout = gzip.open(str(o), 'wt')
            else:
                fout = open(str(o), mode = 'wt')
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

    See :py:mod:`contextlib` for more information.
    """
    if isinstance(files, str):
        raise ValueError(f"This function accepts a list of files as its first \
        argument. You passed a string, {files}. Pass [{files}] instead to pass a list \
        of length 1.")

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

def get_dqs_from_corrected(readsource, correctedsource):
    """
    Get :class:`kbbq.gatk.applybqsr.ModelDQs` from corrected files.

    :return: model DQs
    :rtype: :class:`kbbq.gatk.applybqsr.ModelDQs`
    :param readsource: source of read, original pairs
        (as generated by :func:`generate_reads_from_files`)
    :type readsource: iterable
    :param correctedsource: source of correctedread, correctedoriginal pairs
        (as generated by :func:`generate_reads_from_files`)
    :type correctedsource: iterable
    """
    utils.print_info("Finding errors from corrected files")
    covariates = kbbq.covariate.CovariateData()
    for (read, original), (corrected, original_c) in zip(readsource, correctedsource):
        read.errors = find_corrected_sites(read, corrected)
        covariates.consume_read(read)
    dqs = kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)
    return dqs

def load_subsampled_hash(readsource, ksize, memory, alpha):
    """
    Load the initial subsampled hash.

    :param readsource: source of read, original pairs
        (as generated by :func:`generate_reads_from_files`)
    :type readsource: iterable
    :param ksize: k-mer size to use for hash
    :type ksize: int
    :param memory: memory to use for hash
    :type memory: str
    :param alpha: sampling rate
    :type alpha: double
    :return: subsampled k-mer hash and thresholds
    :rtype: (:class:`khmer.Nodegraph` , :class:`numpy.ndarray` (int))
    """
    utils.print_info("Loading hash")
    graph = kbbq.bloom.create_empty_nodegraph(ksize = ksize, max_mem = memory)

    # 1st pass: load hash
    for read, original in readsource:
        kbbq.bloom.count_read(read, graph, sampling_rate = alpha)

    #calculate and report stats
    fpr = khmer.calc_expected_collisions(graph, force = False, max_false_pos = .15)
    utils.print_info("False positive rate:", str(fpr))
    p_added = kbbq.bloom.p_kmer_added(sampling_rate = alpha, graph = graph)
    utils.print_info("Probability any k-mer was added:", str(p_added))
    thresholds = kbbq.bloom.calculate_thresholds(p_added, graph.ksize())
    utils.print_info("Error thresholds:", thresholds)
    return graph, thresholds

def find_trusted_kmers(readsource, graph, ksize, memory, thresholds):
    """
    Load the trusted k-mer hash.

    :param readsource: source of read, original pairs
        (as generated by :func:`generate_reads_from_files`)
    :type readsource: iterable
    :param graph: The subsampled k-mer graph
    :type graph: :class:`khmer.Nodegraph`
    :param ksize: k-mer size to use for hash
    :type ksize: int
    :param memory: memory to use for hash
    :type memory: str
    :param thresholds: error threshold array
    :type thresholds: :class:`numpy.ndarray` of int
    :return: trusted k-mer hash
    :rtype: :class:`khmer.Nodegraph`
    """
    trustgraph = kbbq.bloom.create_empty_nodegraph(ksize = ksize, max_mem = memory)
    for read, original in readsource:
        errors = kbbq.bloom.infer_read_errors(read, graph, thresholds)
        #don't trust bad quality
        errors[read.qual <= 6] = True
        read.errors = errors
        #find k-size blocks of non-errors and add to the trusted graph
        kbbq.bloom.add_trusted_kmers(read, trustgraph)
    return trustgraph

def fill_read_errors(read, trustgraph):
    """
    Fill the errors attribute of the given :class:`ReadData` object.

    This will modify the input read.

    :param read: read to modify
    :type read: :class:`ReadData`
    :param trustgraph: :class:`khmer.Nodegraph`
    """
    original_seq = read.seq.copy()
    trusted_kmers = kbbq.bloom.kmers_in_graph(read,trustgraph)
    if np.all(~trusted_kmers): #if there are no trusted kmers
        #we fix one if we can and use that as the anchor
        length, base, pos = kbbq.bloom.fix_one(read.seq, trustgraph)
        if length > 0:
            read.seq[pos] = base
            read.errors[pos] = True
        else:
            # if we can't find anything, stop
            covariates.consume_read(read)
            return
    #find errors. read.seq will be altered in the process
    errors, multiple = kbbq.bloom.infer_errors_from_trusted_kmers(read, trustgraph)
    read.errors[errors] = True
    if np.any(read.errors): #if we have errors, check if we have too many
        trusted_sites = np.zeros(len(read), dtype = np.bool)
        t = kbbq.bloom.rolling_window(trusted_sites, trustgraph.ksize())
        t[trusted_kmers,:] = True
        adjust = False
        if not np.any(np.logical_and(trusted_sites, original_seq != read.seq)):
            if not multiple:
                adjust = True
        read.errors = kbbq.bloom.fix_overcorrection(read, ksize, adjust = adjust)
    read.seq = original_seq
    return

def get_dqs_with_hashes(readsource, trustgraph):
    """
    Get DQs from a readsource and trustgraph.

    :param readsource: source of read, original pairs
        (as generated by :func:`generate_reads_from_files`)
    :type readsource: iterable
    :param trustgraph: trusted k-mer hash
    :type trustgraph: :class:`khmer.Nodegraph`
    :return: model DQs
    :rtype: :class:`kbbq.gatk.applybqsr.ModelDQs`
    """
    utils.print_info("Finding errors and building model")
    covariates = kbbq.covariate.CovariateData()
    for read, original in readsource:
        fill_read_errors(read, trustgraph)
        covariates.consume_read(read)
    dqs = kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)
    return dqs

def recalibrate_and_write(read, original, dqs, output, set_oq):
    """
    Recalibrate a read and write it to output given the model parameters.

    :param read: the read to recalibrate
    :type read: :class:`kbbq.read.ReadData`
    :param original: original read data
    :type original: :class:`pysam.AlignedSegment` or :class:`pysam.FastqProxy`
    :param dqs: model parameters
    :type dqs: :class:`kbbq.gatk.applybqsr.ModelDQs`
    :param set_oq: whether OQ tag should be set for SAM records
    :type set_oq: bool
    :param output: file to write to
    :type output: :term:`file object`
    """
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
        output.write(original)
    elif isinstance(original, pysam.FastxRecord):
        original.quality = utils.nparray_to_qualitystring(recalibrated_quals)
        output.write(str(original) + '\n')
    else:
        raise ValueError("Unknown read type {}".format(type(original)))


def recalibrate(files, output, corrected, infer_rg = False, use_oq = False, set_oq = False, ksize = 32, memory = '3G', alpha = .1, gatkreport = None):
    #Check that inputs and outputs match.
    if len(files) != len(output):
        raise ValueError('One output must be specified for each input.')
    validate_files(files)
    bams = load_headers_from_bams(files)

    if corrected != []:
        if len(files) != len(corrected):
            raise ValueError('One corrected file must be specified for each input.')
        validate_files(corrected)
        corrected_bams = [opens_as_bam(i) for i in corrected]

    if gatkreport is None or not pathlib.Path(gatkreport).is_file():
        #gatkreport not provided or the provided report doesn't exist
        #thus we need to train the model and output to file if specified
        if corrected != []:
            #get errors from a corrected file
            with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata, generate_reads_from_files(corrected, corrected_bams, infer_rg, use_oq) as correcteadreaddata: #a list of ReadData generators
                    dqs = get_dqs_from_corrected(itertools.chain.from_iterable(allreaddata), itertools.chain.from_iterable(correcteadreaddata))
        else:
            #get errors with a hash-based approach
            #get_dqs_from_hash()
            #load_subsampled_hash()

            # 1st pass: load hash
            with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData generators
                readsource = itertools.chain.from_iterable(allreaddata) #a single ReadData generator
                graph, thresholds = load_subsampled_hash(readsource, ksize, memory, alpha)

            # 2nd pass: find trusted kmers
            #find_trusted_kmers()
            utils.print_info("Finding trusted k-mers")
            trustgraph = kbbq.bloom.create_empty_nodegraph(ksize = ksize, max_mem = memory)
            with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData generators
                readsource = itertools.chain.from_iterable(allreaddata) #a single ReadData generator
                trustgraph = find_trusted_kmers(readsource, graph, ksize, memory, thresholds)

            # 3rd pass: find errors
            with generate_reads_from_files(files, bams, infer_rg, use_oq) as allreaddata: #a list of ReadData generators
                readsource = itertools.chain.from_iterable(allreaddata) #a single ReadData generator
                dqs = get_dqs_with_hashes(readsource, trustgraph)
        if gatkreport is not None:
            #if gatkreport doesn't exist, save the model to it
            report = kbbq.gatk.bqsr.vectors_to_report(*kbbq.gatk.applybqsr.get_modelvecs_from_covariates(covariates), rg_order = list(kbbq.read.ReadData.rg_to_pu.values()))
            report.write(gatkreport)

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
        #4th pass: recalibrate
        for i, o in zip(allreaddata, opened_outputs):
            for read, original in i:
                recalibrate_and_write(read, original, dqs, o, set_oq)

    utils.print_info("Done")
