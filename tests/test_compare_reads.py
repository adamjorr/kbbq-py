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

def test_bam_to_report(report, recalibratedbam, variable_sites):
    bamreport = compare_reads.bam_to_report(recalibratedbam, 'tests/data/ref.fa', variable_sites)
    assert report.version == bamreport.version
    for s, o in zip(report.tables, bamreport.tables):
        assert s.title == o.title
        assert s.description == o.description
        #assert s.data.equals(o.data) #this is a known issue with floats
        assert_frame_equal(s.data, o.data)

def test_fastq_calibration():
    #TODO: I think the best way to do this is to take the uncalibrated bam reads,
    #instantiate a fastqread object from it, and make sure the resulting calibrated
    #read matches the calibrated bam read.
    pass

# def test_report_creation(report, recalibratedbam, variable_sites):
#     """
#     Test our method to make the recalibrationreport from a bam matches GATK
#     """




