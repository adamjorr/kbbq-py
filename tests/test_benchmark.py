import pytest
import kbbq.main

def test_benchmarks_match(monkeypatch, capfd):
    """
    This is pretty much a sanity check, so we should
    do a nocover and add more thorough tests later.
    """
    import sys
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + "benchmark -b tests/data/conf_regions.bam -r tests/data/ref.fa -v tests/data/conf_regions.vcf.gz --label=label --use-oq".split())
        kbbq.main.main()
    bam_captured = capfd.readouterr()

    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', [sys.argv[0]] + "benchmark -b tests/data/conf_regions.bam -r tests/data/ref.fa -v tests/data/conf_regions.vcf.gz --label=label -f tests/data/reads.fq".split())
        kbbq.main.main()
    fastq_captured = capfd.readouterr()
    assert bam_captured.out == fastq_captured.out