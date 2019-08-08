import pytest
from kbbq import recaltable
import filecmp

###########

#GATKReport

###########

def test_GATKReport_init():
    t = recaltable.GATKReport([])
    assert t.tables == []

def test_GATKReport_fromfile(report_and_file, tmp_path):
    file = report_and_file[1]
    r = recaltable.GATKReport.fromfile(file)
    assert len(r.tables) == 5

    #test error with truncated file
    truncated = tmp_path / 'truncated.txt'
    with open(file, 'r') as f, open(truncated, 'w') as t:
        for i in range(0,320):
            t.write(next(f)) #write lines 0-319
    with pytest.raises(ValueError):
        recaltable.GATKReport.fromfile(truncated)

def test_GATKReport_get_headerstring():
    r = recaltable.GATKReport([])
    assert r.get_headerstring() == '#:GATKReport.v1.1:0'
    r.tables = [0,1,2]
    assert r.get_headerstring() == '#:GATKReport.v1.1:3'

def test_GATKReport_write(report_and_file, tmp_path):
    report, file = report_and_file
    correct = tmp_path / 'report.txt'
    report.write(correct)
    assert filecmp.cmp(file, correct)

def test_GATKReport_str(report_and_file):
    report, file = report_and_file
    reportlines = str(report).splitlines(keepends = True)
    with open(file) as f:
        for fileline, objectline in zip(f, reportlines):
            assert fileline == objectline

def test_GATKReport_repr():
    report = recaltable.GATKReport([])
    assert repr(report) == '#:GATKReport.v1.1:0\n\n'
    report.tables = [0,1,2]
    assert repr(report) == '#:GATKReport.v1.1:3\n0\n1\n2\n'

def test_GATKReport_eq():
    report = recaltable.GATKReport([0,1,2])
    other = recaltable.GATKReport([0,1,2])
    assert report == other

    assert not report == recaltable.GATKReport([])

    assert not report == 3
