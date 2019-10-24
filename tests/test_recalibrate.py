import pytest
from test_compare_reads import FakeRead, bamread_to_fakeread
import kbbq.__main__
from kbbq import recalibrate
from kbbq import compare_reads
import pysam
import numpy as np

def test_find_corrected_sites(simple_fastq_reads):
    for r in simple_fastq_reads:
        r2 = pysam.FastxRecord(name = r.name, sequence = r.sequence, quality = r.quality)
        edited_seq = list(r2.sequence)
        edited_seq[5] = 'C'
        r2.sequence = ''.join(edited_seq)
        correct = np.zeros(len(edited_seq), dtype = np.bool)
        correct[5] = True
        assert np.array_equal(recalibrate.find_corrected_sites(r,r2), correct)

@pytest.fixture()
def uncorr_and_corr_fastq_files(tmp_path):
    uncorr_fastq = tmp_path / 'uncorr.fq'
    corrected_fastq = tmp_path / 'corr.fq'
    with open(uncorr_fastq,'w') as fu, open(corrected_fastq,'w') as fc:
        r = pysam.FastxRecord(
            name = 'foo',
            sequence = 'ATG',
            quality = '((#') #7, 7, 2
        r2 = pysam.FastxRecord(
            name = r.name,
            sequence = 'ACG',
            quality = r.quality)
        fu.write(str(r))
        fc.write(str(r2))
    return str(uncorr_fastq), str(corrected_fastq)

@pytest.fixture()
def uncorr_and_corr_with_rg(tmp_path):
    uncorr_fastq = tmp_path / 'uncorr_withrg.fq'
    corrected_fastq = tmp_path / 'corr_withrg.fq'
    with open(uncorr_fastq,'w') as fu, open(corrected_fastq,'w') as fc:
        r = pysam.FastxRecord(
            name = 'foo/1_RG:Z:bar',
            sequence = 'ATG',
            quality = '((#') #7, 7, 2
        r2 = pysam.FastxRecord(
            name = r.name,
            sequence = 'ACG',
            quality = r.quality)
        fu.write(str(r))
        fc.write(str(r2))
    return str(uncorr_fastq), str(corrected_fastq)

def test_fastq_to_covariate_arrays(uncorr_and_corr_fastq_files, uncorr_and_corr_with_rg):
    correct_pos_errs = np.zeros((1,43,6))
    correct_pos_total = np.zeros((1,43,6))
    correct_pos_errs[0,7,1] = 1
    correct_pos_total[0,7,0] = 1
    correct_pos_total[0,7,1] = 1
    correct_dinuc_errs = np.zeros((1,43,16))
    correct_dinuc_total = np.zeros((1,43,16))
    correct_dinuc_errs[0, 7, compare_reads.Dinucleotide.dinuc_to_int['AT']] = 1
    correct_dinuc_total[0, 7, compare_reads.Dinucleotide.dinuc_to_int['AT']] = 1
    correct_vectors = [np.array([6]), #meanq
        np.array([1]), #rg
        np.array([2]), #rg
        np.array([[0,0,0,0,0,0,0,1] + [0] * 35]), #q
        np.array([[0,0,0,0,0,0,0,2] + [0] * 35]), #q
        correct_pos_errs, #pos
        correct_pos_total, #pos
        correct_dinuc_errs, #dinuc
        correct_dinuc_total] #diunc

    for a,b in zip(correct_vectors, recalibrate.fastq_to_covariate_arrays(
        uncorr_and_corr_fastq_files)):
            assert np.array_equal(a,b)
    for a,b in zip(correct_vectors, recalibrate.fastq_to_covariate_arrays(
        uncorr_and_corr_with_rg, infer_rg = True)):
            assert np.array_equal(a,b)

#this read is used below
correct_read = pysam.FastxRecord(
        name = 'foo',
        sequence = 'ATG',
        quality = '\'\'#') #6, 6, 2

correct_read_with_rg = pysam.FastxRecord(
        name = 'foo/1_RG:Z:bar',
        sequence = 'ATG',
        quality = '\'\'#')

def test_recalibrate_fastq(uncorr_and_corr_fastq_files, uncorr_and_corr_with_rg, capfd):
    recalibrate.recalibrate_fastq(uncorr_and_corr_fastq_files)
    captured = capfd.readouterr()
    assert captured.out == str(correct_read) + '\n'

    #now test with infer_rg = True
    recalibrate.recalibrate_fastq(uncorr_and_corr_with_rg, infer_rg = True)
    captured = capfd.readouterr()
    assert captured.out == str(correct_read_with_rg) + '\n'

    #TODO: we may want test 1000x this read to see a more realistic example

def test_recalibrate_bam():
    with pytest.raises(NotImplementedError):
        recalibrate.recalibrate_bam(None)

def test_recalibrate(uncorr_and_corr_fastq_files, capfd):
    recalibrate.recalibrate(bam = None, fastq = uncorr_and_corr_fastq_files)
    captured = capfd.readouterr()
    assert captured.out == str(correct_read) + '\n'

    with pytest.raises(NotImplementedError):
        recalibrate.recalibrate(fastq = None, bam = 'foo')

    with pytest.raises(NotImplementedError):
        recalibrate.recalibrate(fastq = None, bam = None, gatkreport = 'foo')

    with pytest.raises(ValueError):
        recalibrate.recalibrate(fastq = None, bam = None, gatkreport = None)

def test_recalibrate_main(uncorr_and_corr_fastq_files, monkeypatch, capfd):
    import sys
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["recalibrate",'-f'] + list(uncorr_and_corr_fastq_files) )
        kbbq.__main__.main()
    captured = capfd.readouterr()
    assert captured.out == str(correct_read) + '\n'

    with pytest.raises(NotImplementedError), monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["recalibrate",'-b', 'foo'])
        kbbq.__main__.main()

    with pytest.raises(NotImplementedError), monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + ["recalibrate",'-b', 'foo', '-g', 'bar'])
        kbbq.__main__.main()
