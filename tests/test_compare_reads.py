import pytest
from kbbq import compare_reads
import numpy as np
import filecmp

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

def test_bam_to_covariate_arrays(report, recalibratedbam, variable_sites):
    *vectors = compare_reads.bam_to_covariate_arrays(recalibratedbam, 'tests/data/ref.fa', variable_sites)
    bamreport = vectors_to_report(*vectors, )

# def test_report_creation(report, recalibratedbam, variable_sites):
#     """
#     Test our method to make the recalibrationreport from a bam matches GATK
#     """




