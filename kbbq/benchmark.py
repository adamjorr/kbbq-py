#!/usr/bin/env python3

"""
Utilities for benchmarking calibration methods.
"""

from kbbq import compare_reads
import pysam

def get_ref_dict(reffilename):
    fasta = pysam.FastaFile(reffilename)
    ref = {chrom : np.array(list(fasta.fetch(reference = chrom)), dtype = np.unicode) for chrom in fasta.references}
    return ref

def get_var_sites(vcf):
    vcf = pysam.VariantFile(vcf)
    d = dict()
    for record in vcf:
        d.setdefault(record.chrom, list()).append(int(record.pos)-1)
    return d

def get_full_skips(refdict, var_sites):
    skips = {chrom: np.zeros(len(refdict[chrom]), dtype = np.bool) for chrom in refdict.keys()}
    for chrom in skips.keys():
        variable_positions = var_sites[chrom]
        skips[chrom][variable_positions] = True
    return skips

def get_bam_readname(read):
    """
    Given a bam read, get a canonicalized name.
    The name is the query name + a suffix
    based on whether it is first or second in pair.
    """
    suffix = ("/2" if read.is_read2 else "/1")
    return read.query_name + suffix

def get_fastq_readname(read):
    """
    Given a fastq read, get a canonicalized name.
    This is based on whether the read is first in pair.
    """
    return read.name.split(sep = '_')[0]

def get_error_dict(bamfile, refdict, fullskips):
    """
    Finds errors given an authoritative BAM, ref, and sites to skip.

    The returned dict has readname keys and tuple(np.array, np.array) of bools as values.
    """
    edict = dict()
    for read in bamfile:
        name = get_read_name(read)
        edict[name] = find_read_errors(read, refdict, fullskips)
    return edict

def find_errors_in_fastq(fqreads, edict):
    """
    Finds errors given a readname -> (errors, skips) dict.
    """
    errskips = [edict.get(get_fastq_readname(read)) for read in fqreads]
    return zip(*errskips) #this will return a list of errors and a list of skips

def calculate_q(errors, quals):
    """
    Calculate Actual Q and Predicted Q given a flat
    array of errors and flat array of quality scores.

    The index of the returned array represents the predicted quality score,
    the value represents the actual quality score or the number of bases.
    """
    numtotal = np.bincount(quals.reshape(-1))
    numerrs = np.bincount(quals.reshape(-1), minlength = len(numtotal))
    nonzero = (numtotal != 0)
    p = np.true_divide(numerrs[nonzero], numtotal[nonzero])
    q = compare_reads.p_to_q(p)
    actual_q = np.zeros(len(numtotal), dtype = np.int)
    actual_q[nonzero] = q
    return actual_q, numtotal

def benchmark_fastq(fqfile, bamfile, ref, varfile):
    var_sites = compare_reads.load_positions(varfile)
    ref = get_ref_dict(fafile)
    fullskips = get_full_skips(ref, var_sites)
    edict = get_error_dict(bamfile, ref, fullskips)
    reads = list(pysam.FastxFile(fqfile))
    errors, skips = find_errors_in_fastq(reads, edict)
    quals = [np.array(r.get_quality_array()) for r in reads]
    #turn list of small arrays into 2d arrays
    errors = np.array(errors)
    skips = np.array(skips)
    quals = np.array(quals)
    #get rid of skipped sites (we can't tell if these are errors or not)
    e = errors[~skips]
    quals = quals[~skips]
    return calculate_q(e, quals) #actual_q, ntotal

def benchmark_bam(bamfile, fafile, varfile):
    """
    TODO
    """
    pass

def print_benchmark(actual_q, label, nbases):
    """
    Print benchmark data to a tab separated table.

    Label should be a string and is used for all values
    in the table. No headers are printed. This means
    the benchmark command can be called many times and each
    table can be concatenated together.
    """
    nonzero = (actual_q != 0)
    predicted_q = np.arange(len(actual_q))[nonzero]
    for pq, aq, nb in zip(predicted_q, actual_q, nbases):
        print(pq, aq, label, nb, sep = "\t")

def benchmark(bamfile, fafile, vcffile, fastq = None, label = None):
    """
    Perform the benchmark and print the results to stdout.

    If a fastqfile is not given, the benchmark will be taken from
    the reads in the BAM file. If a fastqfile is given, it should
    have the same read names as the reads in the BAM file. The
    names will then be used to look up where errors occur in the
    reads.
    """
    bam = pysam.AlignmentFile(bamfile, 'r')
    ref = get_ref_dict(fafile)
    var_sites = get_var_sites(vcffile)
    if fastqfile is not None:
        actual_q, nbases = benchmark_fastq(fastqfile, bam, ref, var_sites)
        label = (fastqfile if label is None else label)
    else:
        actual_q, nbases = benchmark_bam(bam, ref, var_sites)
        label = (bamfile if label is None else label)
    print_benchmark(actual_q, label, nbases)
