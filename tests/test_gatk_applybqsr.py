import pytest
from kbbq.gatk import applybqsr
from kbbq import compare_reads as utils
from kbbq import recaltable
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal
import io
import pysam

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

def test_bamread_cycle_covariates(simple_bam_reads):
    assert np.array_equal(applybqsr.bamread_cycle_covariates(simple_bam_reads[0]),
        np.arange(17))
    assert np.array_equal(applybqsr.bamread_cycle_covariates(simple_bam_reads[1]),
        np.flip(-(np.arange(9)+1)))

def test_bamread_dinuc_covariates(simple_bam_reads):
    dinucs = np.array(['TT','TA','AG','GA','AT','TA','AA','AA','AG','GG','GA','AT','TA','AC','CT','TG'], dtype = np.unicode)
    dinuc_ints = utils.Dinucleotide.vecget(dinucs)
    correct = np.concatenate([[-1], dinuc_ints])
    assert np.array_equal(applybqsr.bamread_dinuc_covariates(simple_bam_reads[0], use_oq = False), correct)

    # np.array(['CA','AG','GC','CG','GG','GC','CA','AT'])
    # reversed : np.array(['TA','AC','CG','GG','GC','CG','GA','AC'])
    dinucs = np.array(['AT','TG','GC','CC','CG','GC','CT','TG'])
    dinuc_ints = utils.Dinucleotide.vecget(dinucs)
    correct = np.flip(np.concatenate([[-1], dinuc_ints]))
    assert np.array_equal(applybqsr.bamread_dinuc_covariates(simple_bam_reads[1], use_oq = False), correct)

def test_recalibrate_bamread():
    read = pysam.AlignedSegment()
    read.query_name = 'foo'
    read.query_sequence = 'ATG'
    read.query_qualities = [7,7,2]
    read.set_tag('OQ','((#') #7,7,2
    read.set_tag('RG',0)
    meanq = np.array([10])
    globaldeltaq = np.array([1])
    qscoredeltaq = np.array([[2,2,2,2,2,2,2,2]])
    positiondeltaq = np.zeros((1,8,6))
    positiondeltaq[0,7,:] = 3
    dinucdeltaq = np.zeros([1,8, 16])
    dinucdeltaq[0,7,:] = 5
    assert np.array_equal(applybqsr.recalibrate_bamread(read, meanq,
        globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, np.array([0]), use_oq = False),
            np.array([21,21,2]))
    assert np.array_equal(applybqsr.recalibrate_bamread(read, meanq,
        globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, np.array([0]), use_oq = True),
            np.array([21,21,2]))

def test_get_delta_qs():
    meanq = np.array([10])
    rg_errs = np.array([0])
    rg_total = np.array([1000])
    q_errs = np.array([[0]])
    q_total = np.array([[1000]])
    pos_errs = np.array([[[0]]])
    pos_total = np.array([[[1000]]])
    dinuc_errs = np.array([[[0]]])
    dinuc_total = np.array([[[1000]]])
    rgdeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq = applybqsr.get_delta_qs(
        meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs,
        dinuc_total)
    assert np.array_equal(rgdeltaq, np.array([3]))
    assert np.array_equal(qscoredeltaq,np.array([[2]]))
    assert np.array_equal(positiondeltaq,np.array([[[1]]]))
    assert np.array_equal(dinucdeltaq,np.array([[[1,0]]]))

def test_bam_recalibration(report, recalibratedbam):
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


