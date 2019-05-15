"""
Utilities for recalibrating reads.
"""

from kbbq import compare_reads
from kbbq import recaltable
import pysam
import numpy as np

def recalibrate_fastq(fastq):
    """
    Recalibrate a FASTQ file given a list containing 1) a fastq and 2) a corrected fastq
    """


def recalibrate_bam(bam, use_oq = False, set_oq = False):
    """
    Not yet implemented.
    """
    raise NotImplementedError('Recalibrating a bam is not yet implemented. \
        Try converting your BAM to a FASTQ file with the samtools fastq command. \
        We welcome pull requests if you\'d like to work on this feature.')

def recalibrate(bam, fastq, use_oq = False, set_oq = False, gatkreport = None):
    if gatkreport is not None:
        raise NotImplementedError('GATKreport reading / creation is not yet supported.')
    elif bam is not None:
        recalibrate_bam(bam, use_oq, set_oq)
    elif fastq is not None:
        recalibrate_fastq(fastq)
