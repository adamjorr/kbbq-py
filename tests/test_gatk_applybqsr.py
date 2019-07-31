import pytest
from kbbq.gatk import applybqsr
from kbbq import compare_reads as utils
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal

def test_recalibrate_bamread(report, recalibratedbam):
    """
    Test that our python implementation of GATK calibration matches GATK
    """
    rg_to_pu = utils.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}
    meanq, *vectors = applybqsr.table_to_vectors(report, list(rg_to_pu.values()))
    dqs = applybqsr.get_delta_qs(meanq, *vectors)
    for read in recalibratedbam:
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        recalibrated_quals = applybqsr.recalibrate_bamread(read, meanq, *dqs, rg_to_int, use_oq = True)
        assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)
