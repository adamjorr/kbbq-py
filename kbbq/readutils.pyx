# distutils: language = c++
# distutils: include_dirs = ["kbbq/include/", "kbbq/src/"]
# distutils: sources = kbbq/src/readutils.cc

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from pysam.libchtslib cimport bam1_t
from pysam.libcfaidx cimport kseq_t

cdef extern from "readutils.hh" namespace "readutils":
    cdef cppclass ReadData:
        ReadData() except +
        ReadData(bam1_t* bamrecord) except +
        ReadData(kseq_t* fastqrecord, string rg, int second, string namedelimiter) except +
        string str_qual()
        string canonical_name()
        int get_rg_int()
        string get_pu()
        vector[bool] not_skipped_errors()


