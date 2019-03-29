import pytest
from kbbq import compare_reads

def test_bam_calibration(report, recalibratedbam):
    """
    Test that our python implementation of GATK calibration matches GATK
    """
    id_to_pu = compare_reads.get_id_to_pu(recalibratedbam)
    vectors = compare_reads.table_to_vectors(report, list(id_to_pu.values()), 151)


