# distutils: language = c++

from readutils import ReadData as CReadData
from libcpp.string cimport string
from libcpp.vector cimport vector
from pysam.libchtslib cimport bam1_t, kseq_t
from pysam.licalignedsegment import AlignedSegment
from pysam.libcfaidx import FastqProxy
import numpy as np

cdef class ReadData:
    cdef CReadData creaddata

    def __cinit__(self, AlignedSegment b = None, FastqProxy f = None, string rg = "", int second = 2, string namedelimiter = "_"):
        if b is not None:
            self.creaddata = CReadData(b._delegate)
        elif f is not None:
            self.creaddata = CReadData(f._delegate, rg, second, namedelimiter)
        else:
            raise ValueError("b or f must be specified")

    def str_qual(self):
        return self.creaddata.str_qual()

    def canonical_name(self):
        return self.creaddata.canonical_name()

    def get_rg_int(self):
        return self.creaddata.get_rg_int()

    def get_pu(self):
        return self.creaddata.get_pu()

    def not_skipped_errors(self):
        return self.creaddata.not_skipped_errors()

