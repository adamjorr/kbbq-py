import pytest
from kbbq import recaltable
import filecmp
import numpy as np
import pandas as pd

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
    assert report == recaltable.GATKReport([0,1,2])
    assert not report == recaltable.GATKReport([4,5,6])
    assert not report == recaltable.GATKReport([])
    assert not report == recaltable.GATKReport([], version = '0.0')
    assert not report == 3

######################

# GATKTable

######################

def test_GATKTable_init():
    table = recaltable.GATKTable(title = 'foo', description = 'bar', data = 'baz')
    assert table.title == 'foo'
    assert table.description == 'bar'
    assert table.data == 'baz'

@pytest.fixture
def extablestr():
    return '''#:GATKTable:6:2:%s:%s:%.4f:%.4f:%d:%.2f:;
#:GATKTable:RecalTable0:
ReadGroup                   EventType  EmpiricalQuality  EstimatedQReported  Observations  Errors 
HJCMTCCXX160113.5.AAGGATGT  M                   22.0000             24.3199        210398  1382.00
HK2WYCCXX160124.1.AAGGATGT  M                   22.0000             24.3994        196298  1391.00'''

def test_GATKTable_fromstring(extablestr):
    table = recaltable.GATKTable.fromstring(extablestr)
    assert table.title == 'RecalTable0'
    assert table.description == ''
    assert table.data.shape == (2, 6)

@pytest.fixture
def extable(extablestr):
    return recaltable.GATKTable.fromstring(extablestr)

def test_GATKTable_parse_fmtstring():
    typedict = recaltable.GATKTable.parse_fmtstring(
        header = ['foo','bar','baz','test'],
        fmtstring = '#:GATKTable:0:4:%d:%.4f:%s:%x:;')
    assert typedict == {
        'foo': np.int64,
        'bar': np.float64,
        'baz': str }

def test_GATKTable_get_fmtstring(extable):
   print(extable)
   assert extable.get_fmtstring() == '#:GATKTable:6:2:%s:%s:%.4f:%.4f:%d:%.2f:;'

def test_GATKTable_get_colfmts(extable):
   assert extable.get_colfmts() == ['%s','%s','%.4f','%.4f','%d','%.2f']

def test_GATKTable_get_titlestring():
   table = recaltable.GATKTable(title = 'foo', description = 'bar', data = '')
   assert table.get_titlestring() == '#:GATKTable:foo:bar'

def test_GATKTable_get_datastring(extable, extablestr):
   assert extable.get_datastring() == '\n'.join(extablestr.splitlines()[2:])

def test_GATKTable_get_nrows(extable):
    assert extable.get_nrows() == 2

def test_GATKTable_get_ncols(extable):
    assert extable.get_ncols() == 6

def test_GATKTable_write(extable, extablestr, tmp_path):
    path = tmp_path / 'table.txt'
    with path.open('w') as f:
        extable.write(f)
    assert path.read_text() == extablestr + '\n'

def test_GATKTable_str(extable, extablestr):
    assert str(extable) == extablestr

def test_GATKTable_repr():
    table = recaltable.GATKTable(title = 'foo', description = 'bar',
        data = pd.DataFrame({'spam':['eggs']}))
    assert repr(table) == "#:GATKTable:1:1:%s:;\n#:GATKTable:foo:bar\n   spam\n0  eggs"

def test_GATKTable_eq(extable):
    full = recaltable.GATKTable(title = 'foo', description = 'bar',
        data = pd.DataFrame({'spam':['eggs']}))
    notitle = recaltable.GATKTable(title = '', description = 'bar',
        data = pd.DataFrame({'spam':['eggs']}))
    nodesc = recaltable.GATKTable(title = 'foo', description = '',
        data = pd.DataFrame({'spam':['eggs']}))
    nodata = recaltable.GATKTable(title = 'foo', description = 'bar',
        data = pd.DataFrame())
    #self equal
    for t in (extable, full, notitle, nodesc, nodata):
        assert t == t
    #not equal to others
    assert extable != full
    for t in (notitle, nodesc, nodata):
        assert full != t
        assert t != full
    assert extable != 5
    assert t != 'foo'

######################################

# RecalibrationReport

######################################

def test_RecalibrationReport_init(report_and_file, tmp_path):
    with pytest.raises(ValueError):
        recaltable.RecalibrationReport([])

    r = report_and_file[0]
    assert len(r.tables) == 5
    for table, title in zip(r.tables,
        ['Arguments', 'Quantized', 'RecalTable0', 'RecalTable1', 'RecalTable2']):
            assert table.title == title
            assert table.get_nrows() > 1
            assert table.get_ncols() > 1
    assert r.tables[0].data.index.names == ['Argument']
    assert r.tables[1].data.index.names == ['QualityScore']
    assert r.tables[2].data.index.names == ['ReadGroup']
    assert r.tables[3].data.index.names == ['ReadGroup', 'QualityScore']
    assert r.tables[4].data.index.names == ['ReadGroup', 'QualityScore', 'CovariateName', 'CovariateValue']

def test_RecalibrationReport_str(report_and_file):
    with open(report_and_file[1], 'r') as f:
        assert str(report_and_file[0]) == f.read()
