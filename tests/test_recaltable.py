import pytest
from kbbq import recaltable

def test_report_readwrite(report_and_file):
    report, file = report_and_file
    reportlines = str(report).splitlines(keepends = True)
    with open(file) as f:
        for fileline, objectline in zip(f, reportlines):
            assert fileline == objectline
