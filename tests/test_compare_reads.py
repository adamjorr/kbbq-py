import pytest
from kbbq import compare_reads
import numpy as np

def test_bam_calibration(report, recalibratedbam):
    """
    Test that our python implementation of GATK calibration matches GATK
    """
    rg_to_pu = compare_reads.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}
    meanq, vectors = compare_reads.table_to_vectors(report, list(rg_to_pu.values()))
    dqs = get_delta_qs(meanq, *vectors)
    for read in bam:
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        recalibrated_quals = recalibrate_bamread(read, meanq, *dqs, rg_to_int, compare_reads.Dinucleotide.dinuc_to_int)
        assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)
