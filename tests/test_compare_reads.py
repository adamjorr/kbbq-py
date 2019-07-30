import pytest
from kbbq import compare_reads
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal

def test_bam_calibration(report, recalibratedbam):
    """
    Test that our python implementation of GATK calibration matches GATK
    """
    rg_to_pu = compare_reads.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}
    meanq, *vectors = compare_reads.table_to_vectors(report, list(rg_to_pu.values()))
    dqs = compare_reads.get_delta_qs(meanq, *vectors)
    for read in recalibratedbam:
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        recalibrated_quals = compare_reads.recalibrate_bamread(read, meanq, *dqs, rg_to_int, compare_reads.Dinucleotide.dinuc_to_int)
        assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)

def test_bam_to_report(report, uncalibratedbam, variable_sites):
    # Some day it might pass
    bamreport = compare_reads.bam_to_report(uncalibratedbam, 'tests/data/ref.fa', variable_sites)
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

def test_bam_to_report_fwd(uncalibratedbam, variable_sites):
    # Some day it might pass
    import subprocess
    import pysam
    from kbbq import recaltable

    newbam = 'tests/data/fwdonly.bam'
    newtable = 'tests/data/fwdonly.recal.txt'
    with pysam.AlignmentFile(newbam,'wb', template = uncalibratedbam) as fout:
        for r in uncalibratedbam:
            if not r.is_reverse:
                fout.write(r)
    pysam.index(newbam)
    subprocess.run( ["gatk", 'BaseRecalibrator', '-I'] + \
        [ newbam ] + \
        "-R tests/data/ref.fa --known-sites tests/data/conf_regions.vcf.gz --use-original-qualities".split(' ') + \
        ["-O",newtable])
    report = recaltable.RecalibrationReport.fromfile(newtable)
    bamreport = compare_reads.bam_to_report(pysam.AlignmentFile(newbam), 'tests/data/ref.fa', variable_sites)
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

def test_bam_to_report_rev(uncalibratedbam, variable_sites):
    # Some day it might pass
    import subprocess
    import pysam
    from kbbq import recaltable

    newbam = 'tests/data/revonly.bam'
    newtable = 'tests/data/revonly.recal.txt'
    with pysam.AlignmentFile(newbam,'wb', template = uncalibratedbam) as fout:
        for r in uncalibratedbam:
            if r.is_reverse:
                fout.write(r)
    pysam.index(newbam)
    subprocess.run( ["gatk", 'BaseRecalibrator', '-I'] + \
        [ newbam ] + \
        "-R tests/data/ref.fa --known-sites tests/data/conf_regions.vcf.gz --use-original-qualities".split(' ') + \
        ["-O", newtable])
    report = recaltable.RecalibrationReport.fromfile(newtable)
    bamreport = compare_reads.bam_to_report(pysam.AlignmentFile(newbam), 'tests/data/ref.fa', variable_sites)
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

def test_bam_to_report_fwd_singles(uncalibratedbam, variable_sites, capfd):
    # Some day it might pass
    import subprocess
    import pysam
    from kbbq import recaltable

    newbam = 'tests/data/fwdsingle.bam'
    newtable = 'tests/data/fwdsingle.recal.txt'
    counter = 0 #shortcut ~70 passing reads
    for r in uncalibratedbam:
        if not r.is_reverse:
            if r.get_tag('RG') == 'HK2WY.5':
                counter = counter + 1
                if counter < 120:
                    continue
                with pysam.AlignmentFile(newbam,'wb', template = uncalibratedbam) as fout:
                    fout.write(r)
                pysam.index(newbam)
                with capfd.disabled():
                    #disable GATK printing to STDOUT and STDERR
                    subprocess.run( ["gatk", 'BaseRecalibrator', '-I'] + \
                        [ newbam ] + \
                        "-R tests/data/ref.fa --known-sites tests/data/conf_regions.vcf.gz --use-original-qualities".split(' ') + \
                        ["-O",newtable], capture_output = True)
                report = recaltable.RecalibrationReport.fromfile(newtable)
                bamreport = compare_reads.bam_to_report(pysam.AlignmentFile(newbam), 'tests/data/ref.fa', variable_sites)
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
                        print('read:', str(r))
                        raise

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

def test_fastq_calibration(report, recalibratedbam):
    rg_to_pu = compare_reads.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}
    meanq, *vectors = compare_reads.table_to_vectors(report, list(rg_to_pu.values()))
    dqs = compare_reads.get_delta_qs(meanq, *vectors)
    for read in recalibratedbam:
        fastqread = bamread_to_fakeread(read)
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        if read.is_reverse:
            gatk_calibrated_quals = np.flip(gatk_calibrated_quals)
        rg = rg_to_int[compare_reads.fastq_infer_rg(fastqread)]
        recalibrated_quals = compare_reads.recalibrate_fastq(fastqread, meanq, *dqs, rg = rg, dinuc_to_int = compare_reads.Dinucleotide.dinuc_to_int, secondinpair = compare_reads.fastq_infer_secondinpair(fastqread))
        assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)

