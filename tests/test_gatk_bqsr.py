import pytest
from kbbq.gatk import bqsr
from kbbq.gatk import applybqsr
import numpy as np
from pandas.util.testing import assert_frame_equal
import kbbq.compare_reads as utils
import pysam

def test_bamread_bqsr_cycle(simple_bam_reads):
    assert np.array_equal(bqsr.bamread_bqsr_cycle(simple_bam_reads[0]),
        np.arange(17))
    correct = np.flip(-(np.arange(9)+1))
    assert np.array_equal(bqsr.bamread_bqsr_cycle(simple_bam_reads[1]), correct)
    softclipped = simple_bam_reads[1]
    softclipped.cigartuples = ((4,2),(0,7)) #9M -> 2S7M
    correct[0:2] = 0
    assert np.array_equal(bqsr.bamread_bqsr_cycle(softclipped), correct)

def test_bamread_bqsr_dinuc(simple_bam_reads):
    dinucs = np.array(['TT','TA','AG','GA','AT','TA','AA','AA','AG','GG','GA','AT','TA','AC','CT','TG'], dtype = np.unicode)
    dinuc_ints = utils.Dinucleotide.vecget(dinucs)
    correct = np.concatenate([[-1], dinuc_ints])
    assert np.array_equal(bqsr.bamread_bqsr_dinuc(simple_bam_reads[0], use_oq = False), correct)

    # np.array(['CA','AG','GC','CG','GG','GC','CA','AT'])
    # reversed : np.array(['TA','AC','CG','GG','GC','CG','GA','AC'])
    dinucs = np.array(['AT','TG','GC','CC','CG','GC','CT','TG'])
    dinuc_ints = utils.Dinucleotide.vecget(dinucs)
    correct = np.flip(np.concatenate([[-1], dinuc_ints]))
    assert np.array_equal(bqsr.bamread_bqsr_dinuc(simple_bam_reads[1], use_oq = False), correct)

    softclipped = simple_bam_reads[1]
    softclipped.cigartuples = ((4,2),(0,7)) #9M -> 2S7M
    correct[0:2] = 0
    assert np.array_equal(bqsr.bamread_bqsr_dinuc(softclipped, use_oq = False), correct)

def test_bam_to_bqsr_covariates(tmp_path, simple_fasta, simple_varsites):
    test_bam = tmp_path / 'test.bam'
    header = {'HD' : {'VN' : '1.6', 'SO' : 'coordinate'},
        'SQ' : [{'SN' : 'ref', 'LN' : 45}],
        'RG' : [{'ID': '0', 'PU' : '0'}]}
    with pysam.AlignmentFile(test_bam, 'wb', header = header) as out:
            readstr = 'clipped\t0\tref\t9\t255\t1M9H\t*\t0\t0\tA\t(' #q = 7
            read = pysam.AlignedSegment.fromstring(readstr, out.header)
            read.set_tag('OQ','(') #7,7,2
            read.set_tag('RG','0')
            out.write(read)
    pysam.index(str(test_bam))

    test = bqsr.bam_to_bqsr_covariates(pysam.AlignmentFile(test_bam, 'rb'),
        simple_fasta, simple_varsites)
    correct_pos_errs = np.zeros((1,43,2))
    correct_pos_total = np.zeros((1,43,2))
    correct_pos_total[0,7,0] = 1
    correct_dinuc_errs = np.zeros((1,43,16))
    correct_dinuc_total = np.zeros((1,43,16))
    correct = [np.array([6]), #meanq, float badness
        np.array([0]), #rg
        np.array([1]), #rg
        np.array([[0] * 43]), #q
        np.array([[0] * 7 + [1] + [0] * 35]), #q
        correct_pos_errs, #pos
        correct_pos_total, #pos
        correct_dinuc_errs, #dinuc
        correct_dinuc_total, #dinuc
        ]
    for t, c in zip(test,correct):
        assert np.array_equal(t,c)

def test_bamread_adaptor_boundary(simple_bam_reads):
    assert bqsr.bamread_adaptor_boundary(simple_bam_reads[0]) == 45
    assert bqsr.bamread_adaptor_boundary(simple_bam_reads[1]) == 5
    notlen = simple_bam_reads[0]
    notlen.tlen = 0
    assert bqsr.bamread_adaptor_boundary(notlen) is None
    malformed_rev = simple_bam_reads[1]
    malformed_rev.next_reference_start = 1000
    assert bqsr.bamread_adaptor_boundary(malformed_rev) is None
    malformed_fwd = simple_bam_reads[0]
    malformed_fwd.reference_start = 1000
    assert bqsr.bamread_adaptor_boundary(malformed_fwd) is None

def test_trim_bamread(simple_bam_reads, monkeypatch):
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[0]), np.zeros(17))
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[1]), np.zeros(9))
    notlen = simple_bam_reads[0]
    notlen.tlen = 0
    assert np.array_equal(bqsr.trim_bamread(notlen), np.zeros(17))
    #boundary at beginning on rev read
    monkeypatch.setattr(bqsr, 'bamread_adaptor_boundary', lambda x: 36)
    correct = np.zeros(9, dtype = np.bool)
    correct[0] = True
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[1]), correct)
    #boundary at end of fwd read
    monkeypatch.setattr(bqsr, 'bamread_adaptor_boundary', lambda x: 21)
    correct = np.zeros(17, dtype = np.bool)
    correct[-1] = True
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[0]), correct)
    #boundary left of an insertion
    monkeypatch.setattr(bqsr, 'bamread_adaptor_boundary', lambda x: 13)
    correct = np.zeros(17, dtype = np.bool)
    correct[7:] = True
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[0]), correct)
    #boundary right of an insertion
    monkeypatch.setattr(bqsr, 'bamread_adaptor_boundary', lambda x: 14)
    correct = np.zeros(17, dtype = np.bool)
    correct[10:] = True
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[0]), correct)
    #boundary in deletion
    monkeypatch.setattr(bqsr, 'bamread_adaptor_boundary', lambda x: 18)
    correct = np.zeros(17, dtype = np.bool)
    correct[-3:] = True
    assert np.array_equal(bqsr.trim_bamread(simple_bam_reads[0]), correct)

def test_quantize():
    q_errs = np.array([[1]])
    q_total = np.array([[1000]]) #1/1000 = Q30
    correct = np.array([0] + [93] * 93 )
    assert np.array_equal(bqsr.quantize(q_errs, q_total), correct)

def test_vectors_to_report(report, uncalibratedbam):
    rgs = list(utils.get_rg_to_pu(uncalibratedbam).values())
    correct_vectors = applybqsr.table_to_vectors(report, rg_order = rgs)
    test_report = bqsr.vectors_to_report(*correct_vectors, rg_order = rgs)
    for correct, test in zip(report.tables, test_report.tables):
        if correct.title != 'Quantized':
            assert_frame_equal(correct.data, test.data)

def test_bam_to_report(report, uncalibratedbam, variable_sites):
    bamreport = bqsr.bam_to_report(uncalibratedbam, 'tests/data/ref.fa', variable_sites)
    assert report.version == bamreport.version
    for s, o in zip(report.tables, bamreport.tables):
        assert s.title == o.title
        assert s.description == o.description
        if s.title == 'Quantized':
            continue
        #assert s.data.equals(o.data) #this is a known issue with floats
        try:
            assert_frame_equal(s.data, o.data)
        except AssertionError:
            print('gatk report', s)
            print('bam_to_report', o)
            raise
