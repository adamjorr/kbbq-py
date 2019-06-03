"""
Utilities for benchmarking calibration methods.
"""

from kbbq import compare_reads
import pysam
import numpy as np

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

# def find_errors_in_fastq(fqreads, edict):
#     """
#     Finds errors given a readname -> (errors, skips) dict.
#     """
#     errskips = [edict[get_fastq_readname(read)] for read in fqreads]
#     return zip(*errskips) #this will return a list of errors and a list of skips

def calculate_q(errors, quals):
    """
    Calculate Actual Q and Predicted Q given a flat
    array of errors and flat array of quality scores.

    The index of the returned array represents the predicted quality score,
    the value represents the actual quality score or the number of bases.
    """
    numtotal = np.bincount(quals.reshape(-1))
    numerrs = np.bincount(quals[errors].reshape(-1), minlength = len(numtotal))
    nonzero = (numtotal != 0)
    p = np.true_divide(numerrs[nonzero], numtotal[nonzero])
    q = compare_reads.p_to_q(p)
    actual_q = np.zeros(len(numtotal), dtype = np.int)
    actual_q[nonzero] = q
    return actual_q, numtotal

def benchmark_fastq(fqfile, bamfile, ref, var_sites, bedfh = None):
    fullskips = get_full_skips(ref, var_sites, bedfh)
    edict = get_error_dict(bamfile, ref, fullskips)
    errors, skips, quals = zip(*(edict[get_fastq_readname(r)] + (np.array(r.get_quality_array()),) for r in pysam.FastxFile(fqfile)))
    #turn list of small arrays into 2d arrays
    errors = np.array(errors)
    skips = np.array(skips)
    quals = np.array(quals)
    #get rid of skipped sites (we can't tell if these are errors or not)
    
    e = errors[~skips]
    q = quals[~skips]
    return calculate_q(e, q) #actual_q, ntotal

def get_bamread_quals(read, use_oq = False):
    """
    Return the qualities of the read as an np.array.

    Use the OQ tag for read quals if use_oq = True.
    """
    if use_oq:
        return compare_reads.bamread_get_oq(read)
    else:
        return np.array(read.query_qualities, dtype = np.int)

def benchmark_bam(bamfile, ref, var_sites, use_oq = False, bedfh = None):
    fullskips = get_full_skips(ref, var_sites, bedfh)
    errors, skips, quals = zip(*(compare_reads.find_read_errors(read, ref, fullskips) + (get_bamread_quals(read, use_oq),) for read in bamfile)) #generator
    #turn list of 1d arrays into 2d arrays
    errors = np.array(errors)
    skips = np.array(skips)
    quals = np.array(quals)
    #get rid of skipped sites
    e = errors[~skips]
    q = quals[~skips]
    return calculate_q(e, q)

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

def benchmark(bamfile, fafile, vcffile, fastqfile = None, label = None, use_oq = False, bedfh = None):
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
        actual_q, nbases = benchmark_fastq(fastqfile, bam, ref, var_sites, bedfh)
        label = (fastqfile if label is None else label)
    else:
        actual_q, nbases = benchmark_bam(bam, ref, var_sites, use_oq, bedfh)
        label = (bamfile if label is None else label)
    print_benchmark(actual_q, label, nbases)
