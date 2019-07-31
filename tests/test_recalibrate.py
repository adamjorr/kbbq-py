import pytest
from test_compare_reads import FakeRead, bamread_to_fakeread
from kbbq import compare_reads
from kbbq import recalibrate
from kbbq import benchmark
from kbbq.gatk import bqsr
import pysam
import contextlib
import numpy as np

def test_fastq_to_covariate_arrays(recalibratedbam, variable_sites, monkeypatch):
    """
    TODO
    """
    pass
