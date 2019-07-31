import pytest
from test_compare_reads import FakeRead, bamread_to_fakeread
from kbbq import compare_reads
from kbbq import recalibrate
from kbbq import benchmark
import pysam
import contextlib
import numpy as np

def test_fastq_to_covariate_arrays(recalibratedbam, variable_sites, monkeypatch):
    """
    The bam vectors are definitely right because they create the same recalibration
    table that GATK makes. So if we turn the reads into fake FASTQ reads, the result
    vectors and emulate the errors, the vectors should exactly match.
    """
    reffile = 'tests/data/ref.fa'
    rg_to_pu = compare_reads.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}

    bamvectors = compare_reads.bam_to_covariate_arrays(recalibratedbam, reffile, variable_sites)
    recalibratedbam.reset()
    refdict = benchmark.get_ref_dict(reffile)
    varsites = {chrom : np.array(variable_sites[chrom], dtype = np.int) for chrom in variable_sites.keys()}
    fullskips = benchmark.get_full_skips(refdict, varsites)
    edict = benchmark.get_error_dict(recalibratedbam, refdict, fullskips)
    recalibratedbam.reset()
    bamreads = list(recalibratedbam)
    rgs = {r.get_tag('RG') : None for r in bamreads}
    rg_order = [rg_to_int[r] for r in rgs.keys()]

    @contextlib.contextmanager
    def fakeread_generator(x):
        yield (bamread_to_fakeread(r) for r in bamreads)

    def fake_corrected_sites(uncorr_read, corr_read):
        n = benchmark.get_fastq_readname(uncorr_read)
        e,s = edict[n]
        return e

    with monkeypatch.context() as m:
        m.setattr(pysam, 'FastxFile', fakeread_generator)
        m.setattr(recalibrate, 'find_corrected_sites', fake_corrected_sites) #looks at the read name and gets the errors from a dict
        m.setattr(recalibrate, 'get_fq_skips', lambda x: edict.get(benchmark.get_fastq_readname(x))[1])
        fqvectors = recalibrate.fastq_to_covariate_arrays(('foo','bar'), infer_rg = True)

        # print(bamvectors[5].shape)
        # print(fqvectors[5].shape)
        for b, f in zip(bamvectors, fqvectors):
            print(b[rg_order,...],f)
            assert np.array_equal(b[rg_order,...],f)
