import pytest

@pytest.fixture(params = ['recaltable.txt'])
def recalibration_table(request):
    import kbbq.recaltable
    return kbbq.recaltable.RecalibrationReport.fromfile(request.param)
