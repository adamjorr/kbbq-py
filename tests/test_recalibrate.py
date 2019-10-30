import pytest
from test_compare_reads import FakeRead, bamread_to_fakeread
import kbbq.__main__
import kbbq.read
import kbbq.covariate
from kbbq import recalibrate
from kbbq import compare_reads
import pysam
import numpy as np
import sys
import kbbq.gatk.applybqsr

def test_find_corrected_sites(simple_fastq_reads):
    for fqr in simple_fastq_reads:
        fqr2 = pysam.FastxRecord(name = fqr.name, sequence = fqr.sequence, quality = fqr.quality)
        edited_seq = list(fqr2.sequence)
        edited_seq[5] = 'C'
        fqr2.sequence = ''.join(edited_seq)
        correct = np.zeros(len(edited_seq), dtype = np.bool)
        correct[5] = True
        r = kbbq.read.ReadData.from_fastq(fqr)
        r2 = kbbq.read.ReadData.from_fastq(fqr2)
        assert np.array_equal(recalibrate.find_corrected_sites(r,r2), correct)

@pytest.fixture()
def small_fq_read():
    fqr = pysam.FastxRecord(
        name = 'foo',
        sequence = 'ATG',
        quality = '((#') #7, 7, 2
    return fqr

@pytest.fixture()
def small_fq_file(small_fq_read, tmp_path):
    small_fastq = tmp_path / 'small.fq'
    with small_fastq.open('w') as f:
        f.write(str(small_fq_read))
    return str(small_fastq)

@pytest.fixture()
def small_fq_file_with_rg(small_fq_read, tmp_path):
    small_fastq = tmp_path / 'small_with_rg.fq'
    small_fq_read.name = small_fq_read.name + '/1_RG:Z:bar'
    with small_fastq.open('w') as f:
        f.write(str(small_fq_read))
    return str(small_fastq)

#this read is used below
correct_read = pysam.FastxRecord(
        name = 'foo',
        sequence = 'ATG',
        quality = '\'\'#') #6, 6, 2

correct_read_with_rg = pysam.FastxRecord(
        name = 'foo/1_RG:Z:bar',
        sequence = 'ATG',
        quality = '\'\'#')

def test_recalibrate_fastq(small_fq_file, small_fq_read, tmp_path):
    covariates = kbbq.covariate.CovariateData()
    sfq = kbbq.read.ReadData.from_fastq(small_fq_read, rg = small_fq_file)
    sfq.errors[1] = True
    sfq.skips[sfq.qual < 6] = True
    covariates.consume_read(sfq)
    dqs = kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)
    out = tmp_path / 'out.txt'

    with out.open('w') as f:
        recalibrate.recalibrate_fastq(small_fq_file, dqs, f)
    with out.open('r') as f:
        assert f.read() == str(correct_read) + '\n'

def test_recalibrate_fastq_with_rg(small_fq_file_with_rg, small_fq_read, tmp_path):
    #now test with infer_rg = True
    covariates = kbbq.covariate.CovariateData()
    sfq = kbbq.read.ReadData.from_fastq(small_fq_read)
    sfq.errors[1] = True
    sfq.skips[sfq.qual < 6] = True
    covariates.consume_read(sfq)
    dqs = kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)
    out = tmp_path / 'out.txt'
    with out.open('w') as f:
        recalibrate.recalibrate_fastq(small_fq_file_with_rg, dqs, f, infer_rg = True)
    with out.open('r') as f:
        assert f.read() == str(correct_read_with_rg) + '\n'    

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
