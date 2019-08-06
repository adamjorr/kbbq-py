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

def test_find_corrected_sites(simple_fastq, tmp_path):
    """
    TODO
    """
    corrected_fastq = tmp_path / 'corrected.fq'
    with open(simple_fastq,'r') as fin, open(corrected_fastq,'w') as fout:
        inlines = fin.readlines()
        inlines[5][5] = 'C'
        fout.writelines(inlines)
