import pytest
import kbbq.benchmark as benchmark

def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with \'-m \"not slow\"\')"
    )

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

@pytest.fixture(params = ['tests/data/conf_regions.vcf.gz'])
def variable_sites(request):
    from kbbq import compare_reads
    return compare_reads.get_var_sites(request.param)

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

@pytest.fixture()
def simple_fasta(tmp_path):
    """
    This is given as an example in the SAM spec
    """
    import pysam
    import subprocess
    f = tmp_path / "simple.fa"
    f.write_text(">ref\nAGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT\n")
    # pysam.faidx(str(f))
    subprocess.run(['samtools','faidx',str(f)])
    return str(f)

@pytest.fixture()
def simple_sam(tmp_path):
    """
    This is adapted from the example in the SAM spec
    """
    f = tmp_path / "simple.sam"
    f.write_text("@HD\tVN:1.6\tSO:coordinate\n" +
        "@SQ\tSN:ref\tLN:45\n" +
        "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t==99=?<*+/5:@A99:\n" +
        "r001\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t><>???>>>\tNM:i:1\n")
    return str(f)

@pytest.fixture()
def simple_bam(tmp_path, simple_sam):
    """
    BAM version of simple_sam.
    """
    import pysam
    f = tmp_path / "simple.bam"
    bam = pysam.view("-h","-b",simple_sam)
    f.write_bytes(bam)
    pysam.index(str(f))
    return str(f)

@pytest.fixture()
def simple_vcf(tmp_path):
    """
    Pretends there is a variant at site 10 on the ref.

    This will overlap the read in simple_sam.
    """
    f = tmp_path / "simple.vcf"
    f.write_text('##fileformat=VCFv4.2\n' +
        '##FILTER=<ID=PASS,Description="All filters passed">\n' +
        '##contig=<ID=ref,length=45>\n' +
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' +
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n' +
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsyndip\n" +
        "ref\t10\t.\tG\tT\t30\t.\t.\tGT:AD\t0|1:1,1\n")
    return str(f)

@pytest.fixture()
def simple_bed(tmp_path):
    """
    A simple bed for excluding sites 0-7.

    This should exclude the first 2 sites of r001. 
    """
    f = tmp_path / "simple.bed"
    f.write_text('ref\t8\t46\n')
    return str(f)

@pytest.fixture()
def simple_fastq(tmp_path, simple_bam):
    import pysam
    f = tmp_path / 'simple.fq'
    fq = pysam.fastq('-t','-N','-O',simple_bam)
    fq = fq.replace('\t','_')
    f.write_text(fq)
    return str(f)

@pytest.fixture()
def simple_fastq_reads(simple_fastq):
    import pysam
    return list(pysam.FastxFile(simple_fastq,'r'))

@pytest.fixture()
def simple_bam_reads(simple_bam):
    import pysam
    return list(pysam.AlignmentFile(simple_bam,'rb'))

#################################

# Benchmark fixtures

##################################

@pytest.fixture()
def simple_refdict(simple_fasta):
    return benchmark.get_ref_dict(simple_fasta)

@pytest.fixture()
def simple_varsites(simple_vcf):
    return benchmark.get_var_sites(simple_vcf)

@pytest.fixture()
def simple_bedfh(simple_bed):
    bedfh = open(simple_bed, 'r')
    yield bedfh
    bedfh.close()

@pytest.fixture()
def simple_fullskips(simple_fasta, simple_vcf, simple_bedfh):
    return benchmark.get_full_skips(
        benchmark.get_ref_dict(simple_fasta),
        benchmark.get_var_sites(simple_vcf),
        simple_bedfh)

@pytest.fixture()
def simple_error_dict(simple_bam, simple_refdict, simple_fullskips):
    return benchmark.get_error_dict(pysam.AlignmentFile(simple_bam),
        simple_refdict,
        simple_fullskips)
