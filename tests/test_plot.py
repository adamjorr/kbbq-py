import pytest
import kbbq.plot
import kbbq.main
import sys

def test_benchmark_plot(monkeypatch, benchmarkfile, tmp_path):
    """
    Testing that the output file is actually written to is probably
    the best we can do. No telling if the plot actually makes sense...
    """
    outfile = tmp_path / 'plot.pdf'
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0],'plot', '-o', str(outfile), benchmarkfile])
        kbbq.main.main()
    assert outfile.stat().st_size > 0
