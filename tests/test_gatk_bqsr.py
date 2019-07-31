import pytest
from kbbq.gatk import bqsr
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal

def test_bam_to_report(report, uncalibratedbam, variable_sites):
    bamreport = bqsr.bam_to_report(uncalibratedbam, 'tests/data/ref.fa', variable_sites)
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
