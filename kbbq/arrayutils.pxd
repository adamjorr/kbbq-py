# distutils: include_dirs = include/, src/
# distutils: sources = arrayutils.c

from libcpp cimport bool

cdef extern from "include/arrayutils.h":
    void increment_at_i(int a[], unsigned int i[], unsigned int n)
