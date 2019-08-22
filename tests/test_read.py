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
    assert read.ReadData.rg_to_int[None] == 1
    assert read.ReadData.rg_to_int['foo'] == 2
    assert read.ReadData.numrgs == 3

def test_readdata_from_fastq(simple_fastq_reads):
    fqread = simple_fastq_reads[0]
    r = read.ReadData.from_fastq(fqread, rg = 'foo', second = True)
    assert r.name == 'r001'
    assert read.ReadData.rg_to_int['foo'] == 2
    fqread.name = 'r001/2_RG:Z:foo'
    r = read.ReadData.from_fastq(fqread)
    assert r.name == 'r001'
    assert r.rg == 'foo'
    assert read.ReadData.rg_to_int['foo'] == 2

def test_readdata_load_rgs_from_bamfile():
     bam = pysam.AlignmentFile('-', 'w',
        header = {'RG' : [{'ID': 'FOO', 'SM': 'BAR', 'PU': 'BAZ'}]})
     read.ReadData.load_rgs_from_bamfile(bam)
     assert read.ReadData.rg_to_int['FOO'] == 3
     assert read.ReadData.rg_to_pu['FOO'] == 'BAZ'

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

