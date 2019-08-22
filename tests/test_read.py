import pytest
from kbbq import read
import numpy as np
import pandas as pd
import pysam

# exreaddata = read.ReadData( seq = np.array(['A','T','G']),
#     qual = np.array([6,10,3]),
#     skips = np.array([False, False, True]),
#     name = 'read01',
#     rg = 0,
#     second = False,
#     errors = np.array([False, True, True])
#     )

def test_readdata_init(exreaddata):
    assert read.ReadData.rg_to_pu[0] == 0
    assert read.ReadData.rg_to_int[0] == 0
    assert read.ReadData.numrgs == 1

def test_readdata_from_bamread(simple_bam_reads):
    bamread = simple_bam_reads[0]
    r = read.ReadData.from_bamread(bamread)
    assert np.array_equal(r.qual,
        np.array(pysam.qualitystring_to_array('==99=?<*+/5:@A99:')))
    assert r.rg is None
    bamread.set_tag('OQ', '(' * 17)
    bamread.set_tag('RG','foo')
    r = read.ReadData.from_bamread(bamread, use_oq = True)
    assert np.array_equal(r.qual,np.array([7] * 17))
    assert r.rg == 'foo'
    #the rg 0 from conftest.py has int 0
    assert read.ReadData.rg_to_int[None] == 0
    assert read.ReadData.rg_to_int['foo'] == 1
    assert read.ReadData.numrgs == 2
    bamread.is_reverse = True
    r = read.ReadData.from_bamread(bamread)
    assert np.array_equal(r.qual,
        np.flip(np.array(pysam.qualitystring_to_array('==99=?<*+/5:@A99:'))))
    assert np.array_equal(r.seq,
        np.array(list('CAGTATCCTTTATCTAA')))
    read.ReadData.rg_to_pu = dict()
    read.ReadData.rg_to_int = dict()
    read.ReadData.numrgs = 0

def test_readdata_from_fastq(simple_fastq_reads):
    print(read.ReadData.rg_to_int)
    fqread = simple_fastq_reads[0]
    r = read.ReadData.from_fastq(fqread, rg = 'foo', second = True)
    print(read.ReadData.rg_to_int)
    assert r.name == 'r001'
    assert read.ReadData.rg_to_int['foo'] == 0
    fqread.name = 'r001/2_RG:Z:foo'
    r = read.ReadData.from_fastq(fqread)
    assert r.name == 'r001'
    assert r.rg == 'foo'
    assert read.ReadData.rg_to_int['foo'] == 0
    read.ReadData.rg_to_pu = dict()
    read.ReadData.rg_to_int = dict()
    read.ReadData.numrgs = 0

def test_readdata_load_rgs_from_bamfile():
    bam = pysam.AlignmentFile('-', 'w',
        header = {'RG' : [{'ID': 'FOO', 'SM': 'BAR', 'PU': 'BAZ'}]})
    read.ReadData.load_rgs_from_bamfile(bam)
    assert read.ReadData.rg_to_int['FOO'] == 0
    assert read.ReadData.rg_to_pu['FOO'] == 'BAZ'
    read.ReadData.rg_to_pu = dict()
    read.ReadData.rg_to_int = dict()
    read.ReadData.numrgs = 0

def test_readdata_str_qual(exreaddata):
    assert exreaddata.str_qual() == ["'", "+", "$"]

def test_readdata_canonical_name(exreaddata):
    assert exreaddata.canonical_name() == 'read01/1'

def test_readdata_get_rg_int(exreaddata):
    assert exreaddata.get_rg_int() == 0

def test_readdata_get_pu(exreaddata):
    assert exreaddata.get_pu() == 0

def test_readdata_not_skipped_errors(exreaddata):
    assert np.array_equal(exreaddata.not_skipped_errors(),
        np.array([False, True, False]))

def test_readdata_get_rg_errors(exreaddata):
    errors, valid = exreaddata.get_rg_errors()
    assert np.array_equal(errors, np.array([0], dtype = np.int))
    assert np.array_equal(valid, np.array([0,0], dtype = np.int))

def test_readdata_get_q_errors(exreaddata):
    errors, valid = exreaddata.get_q_errors()
    assert np.array_equal(errors, np.array([10], dtype = np.int))
    assert np.array_equal(valid, np.array([6,10], dtype = np.int))

def test_readdata_get_cycle_array(exreaddata):
    assert np.array_equal(exreaddata.get_cycle_array(), np.array([0,1,2]))
    exreaddata.second = True
    assert np.array_equal(exreaddata.get_cycle_array(), np.array([-1,-2,-3]))

def test_readdata_get_cycle_errors(exreaddata):
    errors, valid = exreaddata.get_cycle_errors()
    assert np.array_equal(errors, np.array([1], dtype = np.int))
    assert np.array_equal(valid, np.array([0,1], dtype = np.int))

def test_readdata_get_dinucleotide_array(exreaddata):
    assert np.array_equal(exreaddata.get_dinucleotide_array(), np.array([-1,1,-1]))

def test_readdata_get_dinuc_errors(exreaddata):
    errors, valid = exreaddata.get_dinuc_errors()
    assert np.array_equal(errors, np.array([1], dtype = np.int))
    assert np.array_equal(valid, np.array([1], dtype = np.int))

def test_readdata_len(exreaddata):
    assert len(exreaddata) == 3

def test_bamread_get_oq(simple_bam_reads):
    bamread = simple_bam_reads[0]
    bamread.set_tag('OQ',''.join(['('] * 17))
    assert np.array_equal(read.bamread_get_oq(bamread), np.array([7] * 17))

def test_bamread_get_quals(simple_bam_reads):
    bamread = simple_bam_reads[0]
    bamread.set_tag('OQ',''.join(['('] * 17))
    assert np.array_equal(read.bamread_get_quals(bamread),
        np.array(pysam.qualitystring_to_array('==99=?<*+/5:@A99:')))
    assert np.array_equal(read.bamread_get_quals(bamread, use_oq = True),
        np.array([7] * 17))
