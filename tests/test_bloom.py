import pytest
import pysam
import numpy as np
import kbbq.bloom

def test_create_empty_nodegraph():
    empty = kbbq.bloom.create_empty_nodegraph(ksize = 3, max_mem = '250M')
    assert not empty.get('ATG')

@pytest.fixture()
def empty_nodegraph():
    return kbbq.bloom.create_empty_nodegraph(ksize = 3, max_mem = '250M')

# exreaddata = read.ReadData( seq = np.array(['A','T','G']),
#     qual = np.array([6,10,3]),
#     skips = np.array([False, False, True]),
#     name = 'read01',
#     rg = 0,
#     second = False,
#     errors = np.array([False, True, True])
#     )

def test_count_read(exreaddata, empty_nodegraph):
    kbbq.bloom.count_read(exreaddata, empty_nodegraph, 1)
    assert empty_nodegraph.get('ATG')


def test_kmers_in_graph(exreaddata, empty_nodegraph):
    kbbq.bloom.count_read(exreaddata, empty_nodegraph, 1)
    assert np.all(kbbq.bloom.overlapping_kmers_in_graph(exreaddata, empty_nodegraph))

def test_overlapping_kmers_in_graph(exreaddata, empty_nodegraph):
    kbbq.bloom.count_read(exreaddata, empty_nodegraph, 1)
    noverlapping = kbbq.bloom.overlapping_kmers_in_graph(exreaddata, empty_nodegraph)
    assert np.all(noverlapping == 1)

def test_overlapping_kmers_possible(exreaddata, empty_nodegraph):
    kbbq.bloom.count_read(exreaddata, empty_nodegraph, 1)
    possible = kbbq.bloom.overlapping_kmers_possible(exreaddata, 3)
    assert np.all(possible == 1)
    
    with pytest.raises(ValueError):
        possible = kbbq.bloom.overlapping_kmers_possible(exreaddata, 4)

def test_p_kmer_added(empty_nodegraph):
    a = kbbq.bloom.p_kmer_added(1, empty_nodegraph)
    assert a == 1
    b = kbbq.bloom.p_kmer_added(.05, empty_nodegraph)
    assert b > 0

def test_calculate_thresholds():
    assert np.all(kbbq.bloom.calculate_thresholds(.05, 32) >= 0)

def test_infer_errors():
    overlapping = np.array([1,2,3,4,5])
    possible = np.array([1,2,3,4,5])
    thresholds = np.array([0,1,2,3,4,4])
    errors = kbbq.bloom.infer_errors(overlapping, possible, thresholds)
    correct = np.array([True, True, True, True, False])
    assert np.array_equal(errors, correct)

def test_infer_read_errors(exreaddata, empty_nodegraph):
    kbbq.bloom.count_read(exreaddata, empty_nodegraph, 1)
    errors = kbbq.bloom.infer_read_errors(exreaddata, empty_nodegraph, np.zeros(100))
    assert np.all(~errors)
    errors = kbbq.bloom.infer_read_errors(exreaddata, empty_nodegraph, np.array([100] * 100))
    assert np.all(errors)

def test_add_trusted_kmers(exreaddata, empty_nodegraph):
    kbbq.bloom.add_trusted_kmers(exreaddata, empty_nodegraph)
    assert not empty_nodegraph.get('ATG')
    exreaddata.errors[:] = False
    kbbq.bloom.add_trusted_kmers(exreaddata, empty_nodegraph)
    assert empty_nodegraph.get('ATG')

def test_find_longest_trusted_block():
    trusted_kmers = np.array([True, True, True])
    a,b = kbbq.bloom.find_longest_trusted_block(trusted_kmers)
    assert a == 0
    assert b == 3
    trusted_kmers = np.array([True, False, True])
    a,b = kbbq.bloom.find_longest_trusted_block(trusted_kmers)
    assert a == 0
    assert b == 1
    trusted_kmers = np.array([False, True, True])
    a,b = kbbq.bloom.find_longest_trusted_block(trusted_kmers)
    assert a == 1
    assert b == 3
    trusted_kmers = np.array([True, True, False])
    a,b = kbbq.bloom.find_longest_trusted_block(trusted_kmers)
    assert a == 0
    assert b == 2

def test_infer_errors_from_trusted_kmers(simple_bam_pairs, empty_nodegraph):
    read = simple_bam_pairs[1][0]
    empty_nodegraph.consume('CAGCGGGAT') #original CAGCGGCAT
    e, m = kbbq.bloom.infer_errors_from_trusted_kmers(read, empty_nodegraph)
    assert e[2] == True

def test_correction_len(simple_bam_pairs, empty_nodegraph):
    read = simple_bam_pairs[1][0]
    empty_nodegraph.consume('CAGCGGGAT') #original CAGCGGCAT
    cor_len, base, m = kbbq.bloom.correction_len(read.seq, empty_nodegraph)
    assert cor_len[0] == 3
    assert base == 'C'
    assert m == False

def test_fix_one(simple_bam_pairs, empty_nodegraph):
    read = simple_bam_pairs[1][0]
    empty_nodegraph.consume('CAGCGGGAT') #original CAGCGGCAT
    cor_len, base, pos = kbbq.bloom.fix_one(read.seq, empty_nodegraph)
    assert cor_len == 6 #this cor len includes ALL the trusted kmers from the beginning.
    assert base == 'C'
    assert pos == 2

# def test_fix_overcorrection():
    # pass

def test_rolling_window():
    a = np.arange(4)
    b = kbbq.bloom.rolling_window(a, 3)
    assert b.shape == (2,3)
    assert np.array_equal(np.array([0,1,2]),b[0,:])
    assert np.array_equal(np.array([1,2,3]),b[1,:])
