cdef extern from "src/arrayutils.c":
    #this block is here to ensure arrayutils gets compiled
    pass

from khmer import Nodegraph
import numpy as np
cimport cython

def increment_at(int[::1,] a, unsigned int[::1,] i):
    #a and i should automatically convert to a memoryview
    increment_at_i(&a[0], &i[0], i.shape[0])
