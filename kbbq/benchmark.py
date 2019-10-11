"""
Utilities for benchmarking calibration methods.
"""

from kbbq import compare_reads
from kbbq import read
from kbbq import recalibrate
import pysam
import numpy as np
import contextlib
import gzip

def get_ref_dict(reffilename):
    fasta = pysam.FastaFile(reffilename)
    ref = {chrom : np.array(list(fasta.fetch(reference = chrom)), dtype = np.unicode) for chrom in fasta.references}
    return ref

def get_var_sites(vcf):
    vcf = pysam.VariantFile(vcf)
    d = dict()
    for record in vcf:
        for i in range(record.start, record.stop, 1):
            d.setdefault(record.chrom, list()).append(i)
    return d

def get_bed_dict(refdict, bedfh):
    beddict = {chrom: np.zeros(len(refdict[chrom]), dtype = np.bool) for chrom in refdict.keys()}
    for bedline in pysam.tabix_iterator(bedfh, parser = pysam.asBed()):
        beddict[bedline.contig][bedline.start:bedline.end] = True #end is 1 past the actual end, so this slice should work properly
    return beddict

def get_full_skips(refdict, var_sites, bedfh = None):
    skips = {chrom: np.zeros(len(refdict[chrom]), dtype = np.bool) for chrom in refdict.keys()}
    for chrom in skips.keys():
        variable_positions = np.array(var_sites[chrom], dtype = np.int)
        skips[chrom][variable_positions] = True

    if bedfh is not None:
        beddict = get_bed_dict(refdict, bedfh)
        for chrom in skips.keys():
            skips[chrom][~beddict[chrom]] = True #skip anything not in the BED

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
        name = get_bam_readname(read)
        e, s = compare_reads.find_read_errors(read, refdict, fullskips)
        #we have to reverse e and s because the samtools fastq command will
        #automatically reverse them, causing the values from the fastq to 
        #not match up.
        if read.is_reverse:
            e = np.flip(e)
            s = np.flip(s)
        edict[name] = (e,s)
    return edict

def calculate_q(numerrs, numtotal):
    """
    Calculate Actual Q and Predicted Q given an array containing counts of errors
    and an array containing total counts for each quality score.

    The index of the returned array represents the predicted quality score,
    the value represents the actual quality score or the number of bases.
    """
    nonzero = (numtotal != 0)
    p = np.true_divide(numerrs[nonzero], numtotal[nonzero])
    q = compare_reads.p_to_q(p)
    actual_q = np.zeros(len(numtotal), dtype = np.int)
    actual_q[nonzero] = q
    return actual_q, numtotal

def benchmark_files(bamfh, ref, fullskips, use_oq = False, fastqfh = None, namedelimiter = '_'):
    """
    Benchmark the given files and return the actual and predicted Q.

    If fastqfh is None (the default), reads from the BAM will be directly benchmarked.
    Otherwise, names from the fastq files will be used to look up
    """
    if fastqfh is not None:
        edict = get_error_dict(bamfile, ref, fullskips)
        readpairs = recalibrate.yield_reads(fastqfh, namedelimiter = namedelimiter)
        error_finder = lambda x, y: edict.get(y.canonical_name()) #x=ReadData, y=read
    else:
        readpairs = recalibrate.yield_reads(bamfh, use_oq = use_oq)
        error_finder = lambda x, y: compare_reads.find_read_errors(x, ref, fullskips)
        
    error_counter = np.zeros(compare_reads.RescaledNormal.maxscore, dtype = np.int64)
    q_counter = np.zeros(compare_reads.RescaledNormal.maxscore, dtype = np.int64)
    for rd, original in readpairs:
        rd.errors, rd.skips = error_finder(rd, original)
        qe, qv = rd.get_q_errors()
        np.add.at(error_counter, qe, 1)
        np.add.at(q_counter, qv, 1)
    return calculate_q(error_counter, q_counter)

def print_benchmark(actual_q, label, nbases):
    """
    Print benchmark data to a tab separated table.

    Label should be a string and is used for all values
    in the table. No headers are printed. This means
    the benchmark command can be called many times and each
    table can be concatenated together.
    """
    nonzero = (nbases != 0)
    nbases = nbases[nonzero]
    predicted_q = np.arange(len(actual_q))[nonzero]
    actual_q = actual_q[nonzero]
    for pq, aq, nb in zip(predicted_q, actual_q, nbases):
        print(pq, aq, label, nb, sep = "\t")

@contextlib.contextmanager
def open_bedfile(bedfile):
    """
    Attempt to open the given bedfile and return a filehandle.

    If bedfile is None, None will be returned.

    This is a context manager; the file will be closed once it goes out of scope.
    """
    if bedfile is None:
        bedfh = None
    else:
        if pathlib.Path(bedfile).suffix == '.gz':
            bedfh = gzip.open(bedfile, 'r')
        else:
            bedfh = open(str(o), mode = 'r')
    try:
        yield bedfh
    finally:
        if bedfh is not None:
            bedfh.close()


def benchmark(bamfile, fafile, vcffile, fastqfile = None, label = None, use_oq = False, bedfile = None, namedelimiter = '_'):
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
    with open_bedfile(bedfile) as bedfh:
        fullskips = get_full_skips(ref, var_sites, bedfh)
    with pysam.AlignmentFile(bamfile, 'r') as bam:
        if fastqfile is not None:
            # bamfh, ref, fullskips, use_oq = False, fastqfh = None, namedelimiter = '_'
            with pysam.FastxFile(fastqfile) as fastqfh:
                actual_q, nbases = benchmark_files(bam, ref, fullskips, use_oq, fastqfh, namedelimiter)
            label = (fastqfile if label is None else label)
        else:
            actual_q, nbases = benchmark_files(bam, ref, fullskips, use_oq)
            label = (bamfile if label is None else label)
    print_benchmark(actual_q, label, nbases)
