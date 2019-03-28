import pytest
from ..kbbq import recaltable

test_report_readwrite(report_and_file):
    report, file = report_and_file
    strfh = io.StringIO()
    with open(file) as f:
        assert f.read() == str(report)
