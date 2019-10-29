import pytest
from kbbq import covariate, read
import numpy as np
import pandas as pd

def test_pad_axis():
    test = np.array([1])
    correct = np.array([1,0,0])
    assert np.array_equal(covariate.pad_axis(test, 0, 2), correct)
    test = np.array([[1]])
    correct = np.array([[1],[0],[0]])
    assert np.array_equal(covariate.pad_axis(test, 0, 2), correct)
    correct = np.array([[1,0,0]])
    assert np.array_equal(covariate.pad_axis(test, 1, 2), correct)

def test_covariate_init():
    test = covariate.Covariate()
    correct = (0,)
    assert test.errors.shape == correct
    assert test.total.shape == correct
    correct = (1,2)
    test = covariate.Covariate(correct)
    assert test.errors.shape == correct
    assert test.total.shape == correct

def test_covariate_pad_axis():
    test = covariate.Covariate()
    test.pad_axis(0)
    assert test.shape() == (1,)
    test = covariate.Covariate((3,4))
    test.pad_axis(1,2)
    assert test.shape() == (3,6)

def test_covariate_pad_axis_to_fit():
    test = covariate.Covariate()
    test.pad_axis_to_fit(0,99)
    assert test.shape() == (100,)
    test = covariate.Covariate((1,2))
    test.pad_axis_to_fit(1,-10)
    assert test.shape() == (1,10)
    test.pad_axis_to_fit(1,0)
    assert test.shape() == (1,10)

def test_covariate_increment():
    test = covariate.Covariate((10,))
    test.increment((0,0),(0,1))
    assert np.array_equal(test.errors, np.zeros([10]))
    assert np.array_equal(test.total, np.array([1] + [0] * 9))

def test_covariate_shape():
    correct = (0,)
    test = covariate.Covariate()
    assert test.shape() == correct
    correct = (1,2)
    test = covariate.Covariate(correct)
    assert test.shape() == correct
    
def test_covariate_getitem():
    test = covariate.Covariate((1,))
    test.increment((0,0),(0,1))
    assert test[0] == (0,1)

def test_covariate_setitem():
    test = covariate.Covariate((1,))
    test[0] = (0,1)
    test.increment((0,0),(0,1))
    assert test[0] == (0,2)

def test_rgcovariate_init():
    test = covariate.RGCovariate()
    assert test.shape() == (0,)

# @pytest.fixture
# def exreaddata():
#     return read.ReadData( seq = np.array(['A','T','G']),
#         qual = np.array([6,10,3]),
#         skips = np.array([False, False, True]),
#         name = 'read01',
#         rg = 0,
#         second = False,
#         errors = np.array([False, True, True])
#         )

def test_rgcovariate_consume_read(exreaddata):
    r = exreaddata

    test = covariate.RGCovariate()
    res = test.consume_read(r)
    rg = read.ReadData.rg_to_int[0]
    assert np.array_equal(res[0], np.array([rg]))
    assert np.array_equal(res[1], np.array([rg] * 2))
    assert test[rg] == (1,2)

def test_rgcovariate_num_rgs():
    test = covariate.RGCovariate()
    assert test.num_rgs() == 0
    test.pad_axis(axis = 0, n = 1)
    assert test.num_rgs() == 1

def test_qcovariate_init():
    test = covariate.QCovariate()
    assert test.shape() == (0,0)

def test_qcovariate_consume_read(exreaddata):
    test = covariate.QCovariate()
    (rge, rgv), (qe, qv) = test.consume_read(exreaddata)
    assert np.array_equal(rge, np.array([0]))
    assert np.array_equal(rgv, np.array([0,0]))
    assert np.array_equal(qe, np.array([10]))
    assert np.array_equal(qv, np.array([6,10]))
    assert test[(0,10)] == (1,1)

def test_qcovariate_num_qs(exreaddata):
    test = covariate.QCovariate()
    assert test.num_qs() == 0
    test.consume_read(exreaddata)
    assert test.num_qs() == 11

def test_cyclecovariate_init():
    test = covariate.CycleCovariate()
    assert test.shape() == (0,0,0)

def test_cyclecovariate_pad_axis():
    test = covariate.CycleCovariate()
    test.pad_axis(axis = 0, n = 1)
    assert test.shape() == (1,0,0)
    with pytest.raises(ValueError):
        test.pad_axis(2,1)
    test.pad_axis(1,1)
    assert test.shape() == (1,1,0)
    test.pad_axis(2,2)
    assert test.shape() == (1,1,2)
    test[(0,0,0)] = (1,1)
    test[(0,0,-1)] = (2,2)
    test.pad_axis(2,4)
    assert test.shape() == (1,1,6)
    print(test.total)
    assert test[0,0,0] == (1,1)
    assert test[0,0,-1] == (2,2)

def test_cyclecovariate_num_cycles():
    test = covariate.CycleCovariate()
    test.pad_axis(2,2)
    assert test.num_cycles() == 1

def test_dinuccovariate_init():
    test = covariate.DinucCovariate()
    assert test.shape() == (0,0,16)

def test_dinuccovariate_num_dinucs():
    test = covariate.DinucCovariate()
    assert test.num_dinucs() == 16

def test_covariatedata_init():
    test = covariate.CovariateData()
    assert test.qcov.shape() == (0,0)
    assert test.cyclecov.shape() == (0,0,0)
    assert test.dinuccov.shape() == (0,0,16)

def test_covariatedata_consume_read(exreaddata):
    test = covariate.CovariateData()
    test.consume_read(exreaddata)
    assert test.qcov.rgcov[0] == (1,2)
    assert test.qcov[0,10] == (1,1)
    assert test.cyclecov[0,6,0] == (0,1)
    assert test.dinuccov[0,10,1] == (1,1)

def test_covariatedata_get_num_rgs():
    test = covariate.CovariateData()
    assert test.get_num_rgs() == 0

def test_covariatedata_get_num_qs():
    test = covariate.CovariateData()
    assert test.get_num_qs() == 0

def test_covariatedata_get_num_cycles():
    test = covariate.CovariateData()
    assert test.get_num_cycles() == 0

def test_covariatedata_get_num_dinucs():
    test = covariate.CovariateData()
    assert test.get_num_dinucs() == 16
