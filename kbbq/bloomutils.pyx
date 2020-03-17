# distutils: language = c++
# distutils: include_dirs = ["/home/adam/code/khmer/include", "/home/adam/code/khmer/third-party/bzip2", "/home/adam/code/khmer/third-party/cqf", "/home/adam/code/khmer/third-party/rollinghash", "/home/adam/code/khmer/third-party/seqan", "/home/adam/code/khmer/third-party/smhasher", "/home/adam/code/khmer/third-party/zlib", "/home/adam/.local/share/virtualenvs/kbbq-py1PUZIt/lib/python3.8/site-packages/numpy/core/include"]

# import khmer._khmer
from khmer._oxli.graphs cimport BoundedCounterType, CpNodegraph, Nodegraph
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
cimport numpy as np
import numpy as np


#python fn
def kmers_in_graph(string s, Nodegraph graph):
    # ingraph = np.zeros(s.length(), dtype = np.bool)
    # cdef bint[::1] ing_view = ingraph
    cdef vector[BoundedCounterType] res
    graph._ng_this.get()[0].get_kmer_counts(s, res)
    # for i in range(s.length() - graph.ksize() - 1):
        # ing_view[i] = res[i]
    ingraph = np.asarray(<BoundedCounterType [:res.size()]>res.data())
    return ingraph

#python fn
def overlapping_kmers_in_graph(string s, Nodegraph graph):
    cdef pair[vector[int], vector[int]] ret = c_overlapping_kmers_in_graph(s, graph._ng_this.get()[0])
    cdef int[:] ingraph_view = <int [:ret.first.size()]>ret.first.data()
    cdef int[:] possible_view = <int [:ret.second.size()]>ret.second.data()
    return np.asarray(ingraph_view), np.asarray(possible_view)

#c++ fn
cdef pair[vector[int], vector[int]] c_overlapping_kmers_in_graph(string s, CpNodegraph graph) nogil:
    """
    Return 2 vectors: (# overlapping kmers, # possible kmers)
    """
    cdef vector[BoundedCounterType] kmers
    graph.get_kmer_counts(s, kmers)
    cdef vector[int] ingraph
    cdef vector[int] possiblevec
    cdef pair[vector[int], vector[int]] ret
    cdef int incount = 0
    cdef int notincount = 0
    cdef int possible = 0
    cdef int ksize = graph.ksize()
    cdef int i


    for i in range(s.length()):
        if i >= ksize:
            # we are past the end of the first kmer
            if kmers[i - ksize]:
                incount = incount - 1
            else:
                notincount = notincount - 1
        if i < kmers.size():
            # we are not past the beginning of the last kmer
            if kmers[i]:
                incount = incount + 1
            else:
                notincount = notincount - 1
        possible = incount + notincount
        ingraph.push_back(incount)
        possiblevec.push_back(possible)

    ret.first = ingraph
    ret.second = possiblevec
    return ret
