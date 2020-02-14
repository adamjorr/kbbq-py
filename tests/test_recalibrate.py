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
def small_fq_read_with_rg(small_fq_read):
    small_fq_read.name = small_fq_read.name + '/1_RG:Z:bar'
    return small_fq_read

@pytest.fixture()
def small_read_with_rg(small_fq_read_with_rg):
    return kbbq.read.ReadData.from_fastq(small_fq_read_with_rg)

@pytest.fixture()
def small_fq_file(small_fq_read, tmp_path):
    small_fastq = tmp_path / 'small.fq'
    with small_fastq.open('w') as f:
        f.write(str(small_fq_read))
    return str(small_fastq)

@pytest.fixture()
def small_read(small_fq_read, small_fq_file):
    return kbbq.read.ReadData.from_fastq(small_fq_read, rg = small_fq_file)

@pytest.fixture()
def small_fq_file_with_rg(small_fq_read_with_rg, tmp_path):
    small_fastq = tmp_path / 'small_with_rg.fq'
    with small_fastq.open('w') as f:
        f.write(str(small_fq_read_with_rg))
    return str(small_fastq)

#this read is used below
@pytest.fixture()
def correct_read():
    return pysam.FastxRecord(
        name = 'foo',
        sequence = 'ATG',
        quality = '\'\'#') #6, 6, 2

@pytest.fixture()
def correct_read_with_rg():
    return pysam.FastxRecord(
        name = 'foo/1_RG:Z:bar',
        sequence = 'ATG',
        quality = '\'\'#')

@pytest.fixture()
def small_fq_dqs(small_read):
    small_read.errors[1] = True
    small_read.skips[small_read.qual < 6] = True
    covariates = kbbq.covariate.CovariateData()
    covariates.consume_read(small_read)
    return kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)

@pytest.fixture()
def small_fq_dqs_with_rg(small_read_with_rg):
    small_read_with_rg.errors[1] = True
    small_read_with_rg.skips[small_read_with_rg.qual < 6] = True
    covariates = kbbq.covariate.CovariateData()
    covariates.consume_read(small_read_with_rg)
    return kbbq.gatk.applybqsr.get_modeldqs_from_covariates(covariates)

def test_recalibrate_read(small_read_with_rg, small_fq_dqs_with_rg, correct_read_with_rg):
    assert np.array_equal(kbbq.recalibrate.recalibrate_read(small_read_with_rg, small_fq_dqs_with_rg),
        np.array(pysam.qualitystring_to_array(correct_read_with_rg.quality)))

def test_recalibrate_fastq(small_fq_file, small_fq_dqs, correct_read, tmp_path):
    out = tmp_path / 'out.txt'
    with out.open('w') as f:
        recalibrate.recalibrate_fastq(small_fq_file, small_fq_dqs, f)
    with out.open('r') as f:
        assert f.read() == str(correct_read) + '\n'

def test_recalibrate_fastq_with_rg(small_fq_file_with_rg, small_fq_dqs_with_rg, correct_read_with_rg, tmp_path):
    #now test with infer_rg = True
    out = tmp_path / 'out.txt'
    with out.open('w') as f:
        recalibrate.recalibrate_fastq(small_fq_file_with_rg, small_fq_dqs_with_rg, f, infer_rg = True)
    with out.open('r') as f:
        assert f.read() == str(correct_read_with_rg) + '\n'    

    #TODO: we may want test 1000x this read to see a more realistic example

def test_find_covariates(small_read):
    correct = kbbq.covariate.CovariateData()
    correct.consume_read(small_read)
    cov = kbbq.recalibrate.find_covariates([[small_read]])
    assert cov == correct

def test_opens_as_bam(small_fq_file, small_report, simple_vcf, simple_bed, simple_bam, simple_sam):
    assert not kbbq.recalibrate.opens_as_bam(small_fq_file)
    assert not kbbq.recalibrate.opens_as_bam(simple_vcf)
    assert not kbbq.recalibrate.opens_as_bam(simple_bed)
    assert kbbq.recalibrate.opens_as_bam(simple_bam)
    assert kbbq.recalibrate.opens_as_bam(simple_sam)

def test_load_headers_from_bams(simple_bam, tmp_path):
    import importlib
    importlib.reload(kbbq.read) #reset read data
    with pysam.AlignmentFile(simple_bam,'rb') as bam: #load bam
        newhead = bam.header.as_dict() #add a read group
        newhead.setdefault('RG',list()).append({'ID':'foo','PU':'bar'})
        with pysam.AlignmentFile(tmp_path/'out.bam','wb', header = newhead) as out:
            for r in bam:
                out.write(r)
    #check rg loaded
    bams = recalibrate.load_headers_from_bams([tmp_path/'out.bam'])
    assert bams[0] == True
    assert kbbq.read.ReadData.rg_to_pu['foo'] == 'bar'
    assert kbbq.read.ReadData.rg_to_int['foo'] == 0
    assert kbbq.read.ReadData.numrgs == 1

def test_validate_files(tmp_path):
    import os
    #raise if does not exist
    shouldraise = tmp_path / "DOESNOTEXIST"
    with pytest.raises(ValueError):
        kbbq.recalibrate.validate_files([shouldraise])
    if os.name == 'posix': #namedpipe only avail on unix
        #raise if exists but isn't a regular file
        namedpipe = tmp_path / "NAMEDPIPE"
        os.mkfifo(namedpipe)
        with pytest.raises(ValueError):
            kbbq.recalibrate.validate_files([namedpipe])
    #should also raise if it's a directory
    with pytest.raises(ValueError):
        kbbq.recalibrate.validate_files([tmp_path])
    #do not raise if exists but is a regular file
    noraise = tmp_path / "exists.txt"
    noraise.write_text("FOO")
    kbbq.recalibrate.validate_files([noraise])

def test_open_outputs(simple_bam, simple_bam_reads, tmp_path):
    files = [simple_bam, tmp_path / "fake.txt", tmp_path / "fake.txt.gz"]
    output = [simple_bam + "out.bam", tmp_path / "fake.out", tmp_path / "fake.out.gz"]
    bams = [True, False, False]
    with kbbq.recalibrate.open_outputs(files, output, bams) as outputs:
        for o in outputs:
            assert not o.closed
        outputs[0].write(simple_bam_reads[0])
        assert outputs[0].header['PG'][-1]['ID'] == 'kbbq'
        outputs[1].write('foo')
        outputs[2].write('bar')
    #check that header adding works properly
    with kbbq.recalibrate.open_outputs([output[0]], [output[0] + "out2.bam"], [True]) as outputs:
        assert outputs[0].header['PG'][-1]['ID'] == 'kbbq.1'

def test_yield_reads(simple_bam, simple_bam_reads, simple_fastq, simple_fastq_reads):
    bamreaddata, bamreads = zip(*list(kbbq.recalibrate.yield_reads(pysam.AlignmentFile(simple_bam))))
    for a,b in zip(bamreads, simple_bam_reads):
        #no __eq__ for AlignedSegment objects :(
        assert str(a) == str(b)

    fqreaddata, fastqreads = zip(*list(kbbq.recalibrate.yield_reads(pysam.FastxFile(simple_fastq))))
    for a,b in zip(fastqreads, simple_fastq_reads):
        #no __eq__ for FastxProxy objects :(
        assert str(a) == str(b)

    with pytest.raises(ValueError):
        foo, bar = kbbq.recalibrate.yield_reads(bamreaddata)

# def test_recalibrate_main(uncorr_and_corr_fastq_files, monkeypatch, capfd):
#     import sys
#     with monkeypatch.context() as m:
#         m.setattr(sys, 'argv', [sys.argv[0]] + ["recalibrate",'-f'] + list(uncorr_and_corr_fastq_files) )
#         kbbq.__main__.main()
#     captured = capfd.readouterr()
#     assert captured.out == str(correct_read) + '\n'

#     with pytest.raises(NotImplementedError), monkeypatch.context() as m:
#         m.setattr(sys, 'argv', [sys.argv[0]] + ["recalibrate",'-b', 'foo'])
#         kbbq.__main__.main()

#     with pytest.raises(NotImplementedError), monkeypatch.context() as m:
#         m.setattr(sys, 'argv', [sys.argv[0]] + ["recalibrate",'-b', 'foo', '-g', 'bar'])
#         kbbq.__main__.main()
