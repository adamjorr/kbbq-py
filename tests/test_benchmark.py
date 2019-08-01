import pytest
import kbbq.main
import kbbq.benchmark as benchmark
import numpy as np
import pysam

def test_get_ref_dict(simple_fasta):
    refdict = {'ref':np.array(list('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT'),dtype = np.unicode)}
    benchdict = benchmark.get_ref_dict(simple_fasta)
    for k in refdict:
        assert k in benchdict
        assert np.array_equal(refdict[k], benchdict[k])

def test_get_var_sites(simple_vcf):
    var_sites = {'ref':[9]}
    assert var_sites == benchmark.get_var_sites(simple_vcf)

def test_get_bed_dict(simple_fasta, simple_bed):
    pos = np.zeros(45, dtype = np.bool)
    pos[8:45] = True
    bed_dict = {'ref': pos}
    with open(simple_bed, 'r') as bedfh:
        benchdict = benchmark.get_bed_dict(benchmark.get_ref_dict(simple_fasta), bedfh)
    for k in bed_dict:
        assert k in benchdict
        assert np.array_equal(bed_dict[k], benchdict[k])

def test_get_full_skips(simple_fasta, simple_vcf, simple_bed):
    skips = np.zeros(45, dtype = np.bool)
    skips[0:8] = True #bed
    skips[9] = True #vcf
    fullskips = {'ref' : skips}
    with open(simple_bed, 'r') as bedfh:
        benchskips = benchmark.get_full_skips(
            benchmark.get_ref_dict(simple_fasta),
            benchmark.get_var_sites(simple_vcf),
            bedfh)
    for k in fullskips:
        assert k in benchskips
        assert np.array_equal(benchskips[k], fullskips[k])

def test_get_bam_readname(simple_bam):
    bam = pysam.AlignmentFile(simple_bam,'rb')
    reads = list(bam)
    assert benchmark.get_bam_readname(reads[0]) == 'r001/1'
    assert benchmark.get_bam_readname(reads[1]) == 'r001/2'

def test_get_fastq_readname(simple_fastq):
    fq = pysam.FastxFile(simple_fastq)
    reads = list(fq)
    assert benchmark.get_fastq_readname(reads[0]) == 'r001/1'
    assert benchmark.get_fastq_readname(reads[1]) == 'r001/2'

def test_get_error_dict(simple_bam, simple_fasta, simple_vcf, simple_bed):
    r1skips = np.zeros(17, dtype = np.bool)
    r1skips[3] = True #from vcf
    r1skips[0:2] = True #from BED
    r2errs = np.zeros(9, dtype = np.bool)
    r2errs[5] = True
    r2errs = np.flip(r2errs) #since r2 is revcomped
    edict = {'r001/1' : (np.zeros(17, dtype = np.bool), r1skips),
        'r001/2' : (r2errs, np.zeros(9, dtype = np.bool))}

    refdict = benchmark.get_ref_dict(simple_fasta)
    with open(simple_bed, 'r') as bedfh:
        fullskips = benchmark.get_full_skips(
            refdict,
            benchmark.get_var_sites(simple_vcf),
            bedfh)

    benchdict = benchmark.get_error_dict(pysam.AlignmentFile(simple_bam),
        refdict,
        fullskips)

    for k in edict:
        assert k in benchdict
        for a,b in zip(benchdict[k], edict[k]):
            assert np.array_equal(a, b)

def test_benchmarks_match(monkeypatch, capfd):
    """
    This is pretty much a sanity check, so we should
    do a nocover and add more thorough tests later.

    It does a FASTQ benchmark and a BAM benchmark and
    checks that they produce the same output, which
    they should since they're the same file.
    """
    import sys
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + "benchmark -b tests/data/conf_regions.bam -r tests/data/ref.fa -v tests/data/conf_regions.vcf.gz --label=label --use-oq".split())
        kbbq.main.main()
    bam_captured = capfd.readouterr()

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + "benchmark -b tests/data/conf_regions.bam -r tests/data/ref.fa -v tests/data/conf_regions.vcf.gz --label=label -f tests/data/allreads.fq".split())
        kbbq.main.main()
    fastq_captured = capfd.readouterr()
    assert bam_captured.out == fastq_captured.out
