import pytest
from kbbq import compare_reads

def test_bam_calibration(report_and_file):
    """
    Test that our python implementation of GATK calibration matches GATK
    """
    vectors = *table_to_vectors(report_and_file[1],)


