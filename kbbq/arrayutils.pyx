cdef extern from "src/arrayutils.c":
    #this block is here to ensure arrayutils gets compiled
    pass

import khmer
import numpy as np
cimport cython

def increment_at(int[::1,] a, unsigned int[::1,] i):
    #a and i should automatically convert to a memoryview
    increment_at_i(&a[0], &i[0], i.shape[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def kmers_in_graph(int[::1,] seq, khmer.Nodegraph graph):
    cdef int ksize = graph.ksize()
    cdef Py_ssize_t nkmers = seq.shape[0] - ksize + 1
    result = np.zeros(nkmers, dtype = np.bool)
    cdef bint[::1,] res_view = result
    for i in range(nkmers):
        result[i] = graph.get_kmer_counts()

def overlapping_kmers_in_graph(int[::1,] seq, bint[::1,] khmer.Nodegraph graph):
    cdef int ksize = graph.ksize()

