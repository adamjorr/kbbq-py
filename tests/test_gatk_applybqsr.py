import pytest
from kbbq.gatk import applybqsr
from kbbq import compare_reads as utils
from kbbq import recaltable
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal
import io

@pytest.fixture
def small_report(report, tmp_path):
    t = """#:GATKReport.v1.1:5
#:GATKTable:2:17:%s:%s:;
#:GATKTable:Arguments:Recalibration argument collection values used in this run
Argument                    Value                                                                   

#:GATKTable:3:94:%d:%d:%d:;
#:GATKTable:Quantized:Quality quantization map
QualityScore  Count    QuantizedScore

#:GATKTable:6:1:%s:%s:%.4f:%.4f:%d:%.2f:;
#:GATKTable:RecalTable0:
ReadGroup  EventType  EmpiricalQuality  EstimatedQReported  Observations  Errors 
1          M                   23.0000              7.0000        200000  1000.00

#:GATKTable:6:1:%s:%d:%s:%.4f:%d:%.2f:;
#:GATKTable:RecalTable1:
ReadGroup  QualityScore  EventType  EmpiricalQuality  Observations  Errors 
1                     7  M                   23.0000        200000  1000.00

#:GATKTable:8:50763:%s:%d:%s:%s:%s:%.4f:%d:%.2f:;
#:GATKTable:RecalTable2:
ReadGroup  QualityScore  CovariateValue  CovariateName  EventType  EmpiricalQuality  Observations  Errors 
1                     7  1               Cycle          M                   23.0000        200000  1000.00
1                     7  AC              Context        M                   23.0000        200000  1000.00

"""
    p = tmp_path / 'small_report.txt'
    p.write_text(t)
    return recaltable.RecalibrationReport.fromfile(p)


def test_table_to_vectors(small_report):
    (meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total,
        dinuc_errs, dinuc_total) = applybqsr.table_to_vectors(small_report, ["1"])
    assert np.array_equal(meanq, np.array([7], dtype = np.float64))
    assert np.array_equal(global_errs, np.array([1000], dtype = np.int64))
    assert np.array_equal(global_total, np.array([200000], dtype = np.int64))
    assert np.array_equal(q_errs, np.array([[0] * 7 + [1000] + [0] * 35], dtype = np.int64))
    assert np.array_equal(q_total, np.array([[0] * 7 + [200000] + [0] * 35], dtype = np.int64))
    correct_pos_errs = np.zeros((1,43,2), dtype = np.int64)
    correct_pos_errs[0,7,0] = 1000
    assert np.array_equal(pos_errs, correct_pos_errs)
    correct_pos_total = np.zeros((1,43,2), dtype = np.int64)
    correct_pos_total[0,7,0] = 200000
    assert np.array_equal(pos_total, correct_pos_total)
    correct_dinuc_errs = np.zeros((1,43,16), dtype = np.int64)
    correct_dinuc_errs[0,7,3] = 1000
    assert np.array_equal(dinuc_errs, correct_dinuc_errs)
    correct_dinuc_total = np.zeros((1,43,16), dtype = np.int64)
    correct_dinuc_total[0,7,3] = 200000
    assert np.array_equal(dinuc_total, correct_dinuc_total)

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
