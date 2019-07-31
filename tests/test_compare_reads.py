import pytest
from kbbq import compare_reads
from kbbq.gatk import applybqsr
import numpy as np
import filecmp
from pandas.util.testing import assert_frame_equal

######################

# FASTQ Recalibration

######################

class FakeRead:
    def __init__(self, name, quality, sequence):
        self.name = name #string
        self.quality = quality #string
        self.sequence = sequence #string

    def get_quality_array(self, offset = 33):
        q = np.array(list(self.quality), dtype = np.unicode)
        quals = np.array(q.view(np.uint32) - offset, dtype = np.uint32)
        return list(quals)

def bamread_to_fakeread(read):
    complement = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    seq = read.query_sequence
    q = read.get_tag('OQ')
    if read.is_reverse:
        seq = ''.join([complement.get(x,'N') for x in reversed(seq)])
        q = q[::-1]
    suffix = ("/2" if read.is_read2 else "/1")
    name = read.query_name + suffix + '_RG:Z:' + read.get_tag('RG')
    return FakeRead(name = name, quality = q, sequence = seq)

def test_fastq_calibration(report, recalibratedbam):
    rg_to_pu = compare_reads.get_rg_to_pu(recalibratedbam)
    rg_to_int = {r:i for i,r in enumerate(rg_to_pu)}
    meanq, *vectors = applybqsr.table_to_vectors(report, list(rg_to_pu.values()))
    dqs = applybqsr.get_delta_qs(meanq, *vectors)
    for read in recalibratedbam:
        fastqread = bamread_to_fakeread(read)
        gatk_calibrated_quals = np.array(read.query_qualities, dtype = np.int)
        if read.is_reverse:
            gatk_calibrated_quals = np.flip(gatk_calibrated_quals)
        rg = rg_to_int[compare_reads.fastq_infer_rg(fastqread)]
        recalibrated_quals = compare_reads.recalibrate_fastq(fastqread, meanq, *dqs, rg = rg, dinuc_to_int = compare_reads.Dinucleotide.dinuc_to_int, secondinpair = compare_reads.fastq_infer_secondinpair(fastqread))
        assert np.array_equal(recalibrated_quals, gatk_calibrated_quals)

