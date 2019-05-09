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

@pytest.fixture(params = ['tests/data/variable_sites.txt'])
def variable_sites(request):
    from kbbq import compare_reads
    return compare_reads.load_positions(request.param)

@pytest.fixture(params = ['tests/data/conf_regions.bam'])
def uncalibratedbam(request):
    import pysam
    return pysam.AlignmentFile(request.param,"r")

@pytest.fixture(params = [
    ["tests/data/conf_regions.bam", "tests/data/ref.fa", "tests/data/conf_regions.vcf.gz"]
    ])
def benchmarkfile(request, monkeypatch, tmp_path, capfd):
    import sys
    import kbbq
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + f"benchmark -b {request.param[0]} -r {request.param[1]} -v {request.param[2]} --label=label --use-oq".split())
        kbbq.main.main()
    p = tmp_path / 'benchmark.txt'
    p.write_text(capfd.readouterr().out)
    return str(p)
