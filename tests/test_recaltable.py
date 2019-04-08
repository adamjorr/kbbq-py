import pytest
from kbbq import recaltable
import filecmp

def test_report_str(report_and_file):
    report, file = report_and_file
    reportlines = str(report).splitlines(keepends = True)
    with open(file) as f:
        for fileline, objectline in zip(f, reportlines):
            assert fileline == objectline

def test_report_readwrite(report_and_file, tmp_path):
    report, file = report_and_file
    cmppath = tmp_path / 'report.txt'
    report.write(str(cmppath))
    assert filecmp.cmp(file, cmppath)
