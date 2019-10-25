import pytest
import kbbq.__main__
import kbbq.benchmark as benchmark
import kbbq.compare_reads
import numpy as np
import pysam

def test_get_ref_dict(simple_refdict):
    correct = {'ref':np.array(list('AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT'),dtype = np.unicode)}
    for k in correct:
        assert k in simple_refdict
        assert np.array_equal(simple_refdict[k], correct[k])

def test_get_var_sites(simple_varsites):
    correct = {'ref':[9]}
    assert correct == simple_varsites

def test_get_bed_dict(simple_fasta, simple_bedfh):
    pos = np.zeros(45, dtype = np.bool)
    pos[8:45] = True
    bed_dict = {'ref': pos}
    benchdict = benchmark.get_bed_dict(benchmark.get_ref_dict(simple_fasta), simple_bedfh)
    for k in bed_dict:
        assert k in benchdict
        assert np.array_equal(bed_dict[k], benchdict[k])

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

def test_get_error_dict(simple_error_dict):
    r1skips = np.zeros(17, dtype = np.bool)
    r1skips[3] = True #from vcf
    r1skips[0:2] = True #from BED
    r2errs = np.zeros(9, dtype = np.bool)
    r2errs[5] = True
    correct = {'r001/1' : (np.zeros(17, dtype = np.bool), r1skips),
        'r001/2' : (r2errs, np.zeros(9, dtype = np.bool))}

    for k in correct:
        assert k in simple_error_dict
        for a,b in zip(simple_error_dict[k], correct[k]):
            assert np.array_equal(a, b)

def test_calculate_q():
    actual = np.array([0, 20, 0, 42], dtype = np.int)
    errors = np.array([0,1,1,0], dtype = np.int)
    total = np.array([0,101,1,1], dtype = np.int)
    a, t = benchmark.calculate_q(errors, total)
    assert np.array_equal(a, actual)
    assert np.array_equal(t, total)

def test_benchmark_files_fastq(simple_fastq, simple_bam, simple_refdict, simple_varsites, simple_fullskips):
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
    # e = np.array([False] * 19 + [True] + [False] * 3, dtype = np.bool)
    e = np.zeros(kbbq.compare_reads.RescaledNormal.maxscore + 1, dtype = np.int)
    e[30] = 1
    t = np.zeros(kbbq.compare_reads.RescaledNormal.maxscore + 1, dtype = np.int)
    np.add.at(t, q, 1)
    # errs = 30
    # Totals: 24 = 3, 28 = 1, 30 = 4, 27 = 2, 9 = 1, 10 = 1
    # 14 = 1, 20 = 1, 25 = 2, 31 = 1, 32 = 1, 29 = 5
    bench_q, bench_t = benchmark.benchmark_files(pysam.AlignmentFile(simple_bam), simple_refdict,
        simple_fullskips, fastqfh = pysam.FastxFile(simple_fastq,'r'))
    correct_q, correct_t = benchmark.calculate_q(e, t)
    print(correct_q)
    print(bench_q)
    assert np.array_equal(bench_q, correct_q)
    assert np.array_equal(bench_t, correct_t)

def test_benchmark_files_bam(simple_bam, simple_refdict, simple_varsites, simple_fullskips):
    #see test_benchmark_files_fastq
    q = np.array([24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] +
        [29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int)
    # e = np.array([False] * 19 + [True] + [False] * 3, dtype = np.bool)
    e = np.zeros(kbbq.compare_reads.RescaledNormal.maxscore + 1, dtype = np.int)
    e[30] = 1
    t = np.zeros(kbbq.compare_reads.RescaledNormal.maxscore + 1, dtype = np.int)
    np.add.at(t, q, 1)
    bench_q, bench_t = benchmark.benchmark_files(pysam.AlignmentFile(simple_bam), simple_refdict,
        simple_fullskips)
    correct_q, correct_t = benchmark.calculate_q(e, t)
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

def test_benchmark(simple_bam, simple_fastq, simple_fasta, simple_vcf, simple_bed, capfd):
    benchmark.benchmark(simple_bam, simple_fasta, simple_vcf, label = 'test',
        bedfile = simple_bed)
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
    benchmark.benchmark(simple_bam, simple_fasta, simple_vcf, fastqfile = simple_fastq,
        label = 'test', use_oq = False, bedfile = simple_bed)
    captured = capfd.readouterr()
    assert captured.out == correct_benchmark

def test_benchmark_main(simple_bam, simple_fasta, simple_vcf, simple_bed, simple_fastq, monkeypatch, capfd):
    import sys
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["benchmark",'-b', simple_bam, '-r', simple_fasta, '-v', simple_vcf, '-d', simple_bed, '--label=test'])
        kbbq.__main__.main()
    captured = capfd.readouterr()
    assert captured.out == correct_benchmark

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["benchmark",'-b', simple_bam, '-r', simple_fasta, '-v', simple_vcf, '-d', simple_bed, '--label=test', '-f', simple_fastq])
        kbbq.__main__.main()
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
        kbbq.__main__.main()
    bam_captured = capfd.readouterr()

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + "benchmark -b tests/data/conf_regions.bam -r tests/data/ref.fa -v tests/data/conf_regions.vcf.gz --label=label -f tests/data/allreads.fq".split())
        kbbq.__main__.main()
    fastq_captured = capfd.readouterr()
    assert bam_captured.out == fastq_captured.out
