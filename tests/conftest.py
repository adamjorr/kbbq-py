import pytest

@pytest.fixture(params = ['tests/data/conf_regions.recal.txt'])
def report_and_file(request):
    from kbbq import recaltable
    return recaltable.RecalibrationReport.fromfile(request.param), request.param

@pytest.fixture
def report(report_and_file):
    return report_and_file[0]

@pytest.fixture(params = ['tests/data/conf_regions.recal.bam'])
def recalibratedbam(request):
    import pysam
    return pysam.AlignmentFile(request.param,"r")
