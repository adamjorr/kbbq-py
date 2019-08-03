import pytest
import kbbq.main
import kbbq.benchmark as benchmark
import numpy as np
import pysam

@pytest.fixture()
def simple_refdict(simple_fasta):
    return benchmark.get_ref_dict(simple_fasta)

def test_get_ref_dict(simple_refdict):
    correct = {'ref':np.array(list('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT'),dtype = np.unicode)}
    for k in correct:
        assert k in simple_refdict
        assert np.array_equal(simple_refdict[k], correct[k])

@pytest.fixture()
def simple_varsites(simple_vcf):
    return benchmark.get_var_sites(simple_vcf)

def test_get_var_sites(simple_varsites):
    correct = {'ref':[9]}
    assert correct == simple_varsites

@pytest.fixture()
def simple_bedfh(simple_bed):
    bedfh = open(simple_bed, 'r')
    yield bedfh
    bedfh.close()

def test_get_bed_dict(simple_fasta, simple_bedfh):
    pos = np.zeros(45, dtype = np.bool)
    pos[8:45] = True
    bed_dict = {'ref': pos}
    benchdict = benchmark.get_bed_dict(benchmark.get_ref_dict(simple_fasta), simple_bedfh)
    for k in bed_dict:
        assert k in benchdict
        assert np.array_equal(bed_dict[k], benchdict[k])

@pytest.fixture()
def simple_fullskips(simple_fasta, simple_vcf, simple_bedfh):
    return benchmark.get_full_skips(
        benchmark.get_ref_dict(simple_fasta),
        benchmark.get_var_sites(simple_vcf),
        simple_bedfh)

def test_get_full_skips(simple_fullskips):
    skips = np.zeros(45, dtype = np.bool)
    skips[0:8] = True #bed
    skips[9] = True #vcf
    fskips = {'ref' : skips}
    for k in fskips:
        assert k in simple_fullskips
        assert np.array_equal(simple_fullskips[k], fskips[k])

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

@pytest.fixture()
def simple_error_dict(simple_bam, simple_refdict, simple_fullskips):
    return benchmark.get_error_dict(pysam.AlignmentFile(simple_bam),
        simple_refdict,
        simple_fullskips)

def test_get_error_dict(simple_error_dict):
    r1skips = np.zeros(17, dtype = np.bool)
    r1skips[3] = True #from vcf
    r1skips[0:2] = True #from BED
    r2errs = np.zeros(9, dtype = np.bool)
    r2errs[5] = True
    r2errs = np.flip(r2errs) #since r2 is revcomped
    correct = {'r001/1' : (np.zeros(17, dtype = np.bool), r1skips),
        'r001/2' : (r2errs, np.zeros(9, dtype = np.bool))}

    for k in correct:
        assert k in simple_error_dict
        for a,b in zip(simple_error_dict[k], correct[k]):
            assert np.array_equal(a, b)

def test_calculate_q():
    errors = np.array([False, True, True] + [False] * 100)
    quals = np.array([3, 2, 1] + [1] * 100, dtype = np.int)
    actual = np.array([0, 20, 0, 42], dtype = np.int)
    total = np.array([0,101,1,1], dtype = np.int)
    a, t = benchmark.calculate_q(errors, quals)
    assert np.array_equal(a, actual)
    assert np.array_equal(t, total)

def test_benchmark_fastq(simple_fastq, simple_bam, simple_refdict, simple_varsites, simple_bedfh):
    #r2errs[5] = True
    # r1skips[3] = True #from vcf
    # r1skips[0:2] = True #from BED
    #r1quals = '==99=?<*+/5:@A99:'
    # [28, 28, 24, 24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25]
    # [24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] excluding skips
    #r2quals = '><>???>>>'
    # [29, 27, 29, 30, 30, 30, 29, 29, 29]
    q = np.array([24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] +
        [29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int)
    e = np.array([False] * 19 + [True] + [False] * 3, dtype = np.bool)
    # errs = 30
    # Totals: 24 = 3, 28 = 1, 30 = 4, 27 = 2, 9 = 1, 10 = 1
    # 14 = 1, 20 = 1, 25 = 2, 31 = 1, 32 = 1, 29 = 5
    bench_q, bench_t = benchmark.benchmark_fastq(simple_fastq, pysam.AlignmentFile(simple_bam), simple_refdict,
        simple_varsites, simple_bedfh)
    correct_q, correct_t = benchmark.calculate_q(e, q)
    assert np.array_equal(bench_q, correct_q)
    assert np.array_equal(bench_t, correct_t)

def test_get_bamread_quals(simple_bam):
    bam = pysam.AlignmentFile(simple_bam,'rb')
    reads = list(bam)
    assert np.array_equal(benchmark.get_bamread_quals(reads[0]),
        np.array([28, 28, 24, 24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25],
            dtype = np.int))
    assert np.array_equal(benchmark.get_bamread_quals(reads[1]),
        np.array([29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int))
    reads[0].set_tag('OQ','!"#$%&\'()*+,-./01')
    reads[1].set_tag('OQ','!"#$%&\'()')
    assert np.array_equal(benchmark.get_bamread_quals(reads[0], use_oq = True), np.arange(17))
    assert np.array_equal(benchmark.get_bamread_quals(reads[1], use_oq = True), np.arange(9))

def test_benchmark_bam(simple_bam, simple_refdict, simple_varsites, simple_bedfh):
    #see test_benchmark_fastq
    q = np.array([24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] +
        [29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int)
    e = np.array([False] * 19 + [True] + [False] * 3, dtype = np.bool)
    bench_q, bench_t = benchmark.benchmark_bam(pysam.AlignmentFile(simple_bam), simple_refdict,
        simple_varsites, bedfh = simple_bedfh)
    correct_q, correct_t = benchmark.calculate_q(e, q)
    assert np.array_equal(bench_q, correct_q)
    assert np.array_equal(bench_t, correct_t)

def test_print_benchmark(capfd):
    actual = np.array([0, 20, 0, 42], dtype = np.int)
    total = np.array([0,101,1,1], dtype = np.int)
    benchmark.print_benchmark(actual, 'test', total)
    captured = capfd.readouterr()
    assert captured.out == "1\t20\ttest\t101\n" +\
        "2\t0\ttest\t1\n" +\
        "3\t42\ttest\t1\n"

correct_benchmark = "9\t42\ttest\t1\n" +\
    "10\t42\ttest\t1\n" +\
    "14\t42\ttest\t1\n" +\
    "20\t42\ttest\t1\n" +\
    "24\t42\ttest\t3\n" +\
    "25\t42\ttest\t2\n" +\
    "27\t42\ttest\t2\n" +\
    "28\t42\ttest\t1\n" +\
    "29\t42\ttest\t5\n" +\
    "30\t6\ttest\t4\n" +\
    "31\t42\ttest\t1\n" +\
    "32\t42\ttest\t1\n"

def test_benchmark(simple_bam, simple_fastq, simple_fasta, simple_vcf, simple_bedfh, capfd):
    benchmark.benchmark(simple_bam, simple_fasta, simple_vcf, label = 'test',
        bedfh = simple_bedfh)
    captured = capfd.readouterr()
    assert captured.out == correct_benchmark
    #q = np.array([24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] +
    #    [29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int)
    #e = np.array([False] * 19 + [True] + [False] * 3, dtype = np.bool)
    #correct_q, correct_t = benchmark.calculate_q(e, q)
    #correct q
    #[ 0,  0,  0,  0,  0,  0,  0,  0,  0, 42, 42,  0,  0,  0, 42,  0,  0,
    #  0,  0,  0, 42,  0,  0,  0, 42, 42,  0, 42, 42, 42,  6, 42, 42]
    #correct t
    #[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
    # 0, 0, 3, 2, 0, 2, 1, 5, 4, 1, 1]
    simple_bedfh.seek(0)
    benchmark.benchmark(simple_bam, simple_fasta, simple_vcf, fastqfile = simple_fastq,
        label = 'test', use_oq = False, bedfh = simple_bedfh)
    captured = capfd.readouterr()
    assert captured.out == correct_benchmark

def test_benchmark_main(simple_bam, simple_fasta, simple_vcf, simple_bed, simple_fastq, monkeypatch, capfd):
    import sys
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["benchmark",'-b', simple_bam, '-r', simple_fasta, '-v', simple_vcf, '-d', simple_bed, '--label=test'])
        kbbq.main.main()
    captured = capfd.readouterr()
    assert captured.out == correct_benchmark

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["benchmark",'-b', simple_bam, '-r', simple_fasta, '-v', simple_vcf, '-d', simple_bed, '--label=test', '-f', simple_fastq])
        kbbq.main.main()
    captured = capfd.readouterr()
    assert captured.out == correct_benchmark

@pytest.mark.slow
def test_benchmarks_main_big(monkeypatch, capfd):
    """
    This does a FASTQ benchmark and a BAM benchmark and
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
