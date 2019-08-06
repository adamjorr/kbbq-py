import pytest
from kbbq import compare_reads
from kbbq.gatk import applybqsr
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal
import pysam

######################

# FASTQ Recalibration

######################

class FakeRead:
    def __init__(self, name, quality, sequence):
        self.name = name #string
        self.quality = quality #string
        self.sequence = sequence #string

    def get_quality_array(self, offset = 33):
        q = np.array(list(self.quality), dtype = np.unicode)
        quals = np.array(q.view(np.uint32) - offset, dtype = np.uint32)
        return list(quals)

def bamread_to_fakeread(read):
    complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    seq = read.query_sequence
    q = read.get_tag('OQ')
    if read.is_reverse:
        seq = ''.join([complement.get(x,'N') for x in reversed(seq)])
        q = q[::-1]
    suffix = ("/2" if read.is_read2 else "/1")
    name = read.query_name + suffix + '_RG:Z:' + read.get_tag('RG')
    return FakeRead(name = name, quality = q, sequence = seq)

@pytest.mark.slow
def test_fastq_calibration(report, recalibratedbam):
    rg_to_pu = compare_reads.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}
    meanq, *vectors = applybqsr.table_to_vectors(report, list(rg_to_pu.values()))
    dqs = applybqsr.get_delta_qs(meanq, *vectors)
    for read in recalibratedbam:
        fastqread = bamread_to_fakeread(read)
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        if read.is_reverse:
            gatk_calibrated_quals = np.flip(gatk_calibrated_quals)
        rg = rg_to_int[compare_reads.fastq_infer_rg(fastqread)]
        recalibrated_quals = compare_reads.recalibrate_fastq(fastqread, meanq, *dqs, rg = rg, dinuc_to_int = compare_reads.Dinucleotide.dinuc_to_int, secondinpair = compare_reads.fastq_infer_secondinpair(fastqread))
        assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)

def test_tstamp():
    import datetime
    correct = datetime.datetime.today()
    correct = correct.replace(microsecond = 0)
    test = datetime.datetime.fromisoformat(compare_reads.tstamp()[2:-2])
    assert correct == test

def test_load_positions(simple_bed):
    correct = {'ref':list(range(8,46))}
    assert compare_reads.load_positions(simple_bed) == correct

def test_get_var_sites(simple_vcf):
    correct = {'ref':[9]}
    assert compare_reads.get_var_sites(simple_vcf) == correct

def test_train_regression():
    q = np.array([28, 28, 24, 24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] #0:16
        + [29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int) #17:
    e = np.zeros(len(q))
    e[21] = True
    lr = compare_reads.train_regression(q, e)
    q = q[:,np.newaxis]
    predicted = lr.predict_proba(q)
    assert predicted.shape[0] == q.shape[0]

def test_regression_recalibrate():
    q = np.array([28, 28, 24, 24, 28, 30, 27,  9, 10, 14, 20, 25, 31, 32, 24, 24, 25] #0:16
        + [29, 27, 29, 30, 30, 30, 29, 29, 29], dtype = np.int) #17:
    e = np.zeros(len(q))
    e[21] = True
    lr = compare_reads.train_regression(q,e)
    newq = compare_reads.regression_recalibrate(lr, q)
    assert newq.shape == q.shape
    assert not np.all(newq == q)

def test_find_read_errors(simple_bam, simple_refdict, simple_fullskips, monkeypatch):
    """
    This will need more testing for some edge cases probably
    """
    import collections
    r1skips = np.zeros(17, dtype = np.bool)
    r1skips[3] = True #from vcf
    r1skips[0:2] = True #from BED
    r2errs = np.zeros(9, dtype = np.bool)
    r2errs[5] = True

    bam = pysam.AlignmentFile(simple_bam,'rb')
    reads = list(bam)
    e, s = compare_reads.find_read_errors(reads[0], simple_refdict, simple_fullskips)
    assert np.array_equal(e,np.zeros(17, dtype = np.bool))
    assert np.array_equal(s,r1skips)
    e, s = compare_reads.find_read_errors(reads[1], simple_refdict, simple_fullskips)
    assert np.array_equal(e,r2errs)
    assert np.array_equal(s,np.zeros(9, dtype = np.bool))

    #test hard clip block
    readstr = 'clipped\t0\tref\t9\t255\t1M9H\t*\t0\t0\tA\t)'
    clippedread = pysam.AlignedSegment.fromstring(readstr,bam.header)
    e, s = compare_reads.find_read_errors(clippedread, simple_refdict, simple_fullskips)
    assert np.array_equal(e, np.array([False]))
    assert np.array_equal(s, np.array([False], dtype = np.bool))

    # test exception for invalid CIGAR
    # though the amount of effort required to do this means it's unlikely
    # to ever happen to a user...
    FakeRead = collections.namedtuple('FakeRead', clippedread.__dir__(), rename = True)
    shallowvalues = {k:getattr(clippedread,k) for k in clippedread.__dir__()}
    shallowvalues['cigartuples'] =  [('L',9)]
    clippedread = FakeRead(*shallowvalues)
    with pytest.raises(ValueError):
        e, s = compare_reads.find_read_errors(clippedread, simple_refdict, simple_fullskips)

def test_RescaledNormal_prior():
    assert compare_reads.RescaledNormal.prior(0) == np.log(.9)
    assert np.array_equal(compare_reads.RescaledNormal.prior(np.arange(43)),
        compare_reads.RescaledNormal.prior_dist)
    assert np.all(compare_reads.RescaledNormal.prior_dist < 0) # 0 <= p <= 1

def test_Dinucleotide_vecget():
    dinucs = ['AA','AT','AG','AC'] + \
        ['TA','TT','TG','TC'] + \
        ['GA','GT','GG','GC'] + \
        ['CA','CT','CG','CC']
    assert dinucs == compare_reads.Dinucleotide.dinucs
    dinucs = np.array(dinucs, dtype = np.unicode)
    assert np.array_equal(compare_reads.Dinucleotide.vecget(dinucs),
        np.arange(16))
    assert compare_reads.Dinucleotide.vecget(np.array(['AA'], dtype = np.unicode)) == 0

def test_gatk_delta_q():
    prior_q = np.array([10,20,30])
    numerrs = np.array([10, 200, 0])
    numtotal = np.array([1000, 1000, 50000])
    dq = compare_reads.gatk_delta_q(prior_q, numerrs, numtotal)
    print(dq)
    assert dq.shape == prior_q.shape
    assert dq[0] > 0
    assert dq[1] < 0
    assert dq[2] > 0
    assert np.all(dq + prior_q <= 42)
    assert np.all(dq + prior_q > 0)

def test_p_to_q():
    p = np.array([.2, .3, .4, .1, .01, .001])
    correct = np.array([6, 5, 3, 10, 20, 30])
    assert np.array_equal(compare_reads.p_to_q(p), correct)
    allq = np.arange(43)
    allp = compare_reads.q_to_p(allq)
    #because of floats we just need to be within 1 of the correct answer
    difference = allq - compare_reads.p_to_q(allp)
    assert np.all(np.logical_and(difference >= 0, difference <= 1)) #should be non negative

def test_q_to_p():
    q = np.array([6, 10, 20, 30])
    correct = np.array([.251188643, .1, .01, .001])
    assert np.allclose(compare_reads.q_to_p(q), correct)
    
def test_generic_cycle_covariate():
    assert np.array_equal(compare_reads.generic_cycle_covariate(17), np.arange(17))
    assert np.array_equal(compare_reads.generic_cycle_covariate(17, True), -(np.arange(17)+1))

def test_generic_dinuc_covariate():
    s = np.array(list('ATGCATGC'))
    q = np.array([10] * 8)
    dinucs = np.array(['AT', 'TG', 'GC', 'CA', 'AT', 'TG', 'GC'], dtype = np.unicode)
    dinuc_ints = compare_reads.Dinucleotide.vecget(dinucs)
    correct = np.concatenate([[-1], dinuc_ints])
    assert np.array_equal(compare_reads.generic_dinuc_covariate(s,q), correct)

    #test N's
    s[1] = 'N'
    correct[1] = -1
    correct[2] = -1
    assert np.array_equal(compare_reads.generic_dinuc_covariate(s,q), correct)

    #test N's + one below minscore
    q[6] = 2
    correct[6] = -1
    assert np.array_equal(compare_reads.generic_dinuc_covariate(s,q), correct)

def test_fastq_cycle_covariates(simple_fastq_reads):
    assert np.array_equal(compare_reads.fastq_cycle_covariates(simple_fastq_reads[0]),
        np.arange(17))
    assert np.array_equal(compare_reads.fastq_cycle_covariates(simple_fastq_reads[1], secondinpair = True),
        -(np.arange(9)+1))

def test_fastq_dinuc_covariates(simple_fastq_reads):
    dinucs = np.array(['TT','TA','AG','GA','AT','TA','AA','AA','AG','GG','GA','AT','TA','AC','CT','TG'], dtype = np.unicode)
    dinuc_ints = compare_reads.Dinucleotide.vecget(dinucs)
    correct = np.concatenate([[-1], dinuc_ints])
    assert np.array_equal(compare_reads.fastq_dinuc_covariates(simple_fastq_reads[0]), correct)

    # np.array(['CA','AG','GC','CG','GG','GC','CA','AT'])
    # reversed : np.array(['TA','AC','CG','GG','GC','CG','GA','AC'])
    dinucs = np.array(['AT','TG','GC','CC','CG','GC','CT','TG'])
    dinuc_ints = compare_reads.Dinucleotide.vecget(dinucs)
    correct = np.concatenate([[-1], dinuc_ints])
    assert np.array_equal(compare_reads.fastq_dinuc_covariates(simple_fastq_reads[1]), correct)

def test_fastq_infer_secondinpair(simple_fastq_reads):
    assert compare_reads.fastq_infer_secondinpair(simple_fastq_reads[0]) == False
    assert compare_reads.fastq_infer_secondinpair(simple_fastq_reads[1]) == True

def test_fastq_infer_rg(simple_fastq_reads):
    for r in simple_fastq_reads:
        r.name = r.name + '_RG:Z:FOO'
        assert compare_reads.fastq_infer_rg(r) == 'FOO'

def test_recalibrate_fastq(simple_fastq_reads):
    read = pysam.FastxRecord(
        name = 'foo',
        sequence = 'ATG',
        quality = '((#') #7, 7, 2
    meanq = np.array([10])
    globaldeltaq = np.array([1])
    qscoredeltaq = np.array([[2,2,2,2,2,2,2,2]])
    positiondeltaq = np.zeros((1,8,6))
    positiondeltaq[0,7,:] = 3
    dinucdeltaq = np.zeros([1,8, 16])
    dinucdeltaq[0,7,:] = 5
    assert np.array_equal(compare_reads.recalibrate_fastq(read, meanq,
        globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, np.array([0]),
        compare_reads.Dinucleotide.dinuc_to_int), np.array([21,21,2]))

def test_bamread_get_oq(simple_bam_reads):
    simple_bam_reads[0].set_tag('OQ','!"#$%&\'()*+,-./01')
    simple_bam_reads[1].set_tag('OQ','!"#$%&\'()')
    assert np.array_equal(compare_reads.bamread_get_oq(simple_bam_reads[0]), np.arange(17))
    assert np.array_equal(compare_reads.bamread_get_oq(simple_bam_reads[1]), np.arange(9))

def test_get_rg_to_pu():
    bam = pysam.AlignmentFile('-', 'w', header = {'RG' : [{'ID': 'FOO', 'SM': 'BAR', 'PU': 'BAZ'}]})
    assert compare_reads.get_rg_to_pu(bam) == {'FOO': 'BAZ'}
