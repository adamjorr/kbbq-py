"""
A data class and helper functions for read covariates.
"""

import abc
from . import read

def pad_start(array, n = 1):
    """
    Pad the beginning axis of :class:`np.ndarray` array with n zeros.
    """
    return np.append(array, [np.zeros(array.shape[1:], dtype = np.int)] * n, axis = 0)

def pad_end(array, n = 1):
    """
    Pad the end axis of :class:`np.ndarray` array with n zeros.
    """
    return np.append(array, np.zeros(array.shape[0:-1] + (n,), dtype = np.int), axis = -1)

class Covariate(metaclass = abc.ABCMeta):
    """
    A base class for covariates.

    Holds 2 :class:`np.ndarray`, :attr:`total` for
    number of observations and :attr:`errors` for number of errors.

    :meth:`__getitem__` will pass the key to both arrays and return a tuple
    of (errors, total). As a pneumonic, you can remember this as a fraction
    from left to right, as total should always be larger than errors.

    :meth:`__setitem__` will also pass the key to both arrays, and will expect
    value to be a tuple of values, setting the first value to errors and the
    second to total.
    """
    def __init__(self, shape = 0):
        self.errors = np.zeros(shape, dtype = np.int)
        self.total = np.zeros(shape, dtype = np.int)
        
    def pad_start(self, n = 1):
        """
        Expand the first dimension by n.
        """
        self.errors = pad_start(self.errors, n)
        self.total = pad_start(self.total, n)
        
    def pad_end(self, n = 1):
        """
        Expand the last dimension by n.
        """
        self.errors = pad_end(self.errors, n)
        self.total = pad_end(self.total, n)

    def pad_to_fit(self, idx):
        """
        Pad the arrays to fit indices idx. Idx should be a tuple corresponding to (start, end)

        Use a 0 value if you don't want to pad the start or end.
        """

    def increment(self, idx, value = (1,1)):
        """
        Increment the given indices by value. idx should be a tuple of values interpretable by :func:`np.add`.

        The first indices given will increment the error array.
        The second indices given will increment the total array.
        """
        np.add.at(self.errors, idx[0], value[0])
        np.add.at(self.total, idx[1], value[1])

    def shape(self):
        """
        Return the shape of the underlying arrays.
        """
        return self.total.shape

    @abc.abstractmethod
    def consume_read(self, read):
        """
        Classes that inherit from this abstract class must implement consume_read.

        This function should accept a :class:`ReadData` and
        increment the arrays properly.
        """
        pass

    def __getitem__(self, key):
        """
        Pass key to each :class:`np.ndarray` and return the result as a tuple.

        The result will be in order of (errors, total). As a pneumonic,
        you can remember this as a fraction from left to right,
        as total should always be larger than errors.
        """
        return (self.errors[key], self.total[key])

    def __setitem__(self, key, value):
        """
        Pass key to each :class:`np.ndarray` and assign each value.

        Value must be a tuple of values to assign. The first value
        is assigned to errors and the second to total.
        """
        self.errors[key] = value[0]
        self.total[key] = value[1]

class RGCovariate(Covariate):
    """
    A 1d covariate array in the read group dimension.
    """
    def __init__(self):
        super().__init__(shape = 0)

    def consume_read(self, read):
        rg = read.get_rg_int()

        #check boundaries
        largest_idx = self.shape()[0] - 1
        if rg > largest_idx:
            padsize = abs(largest_idx - rg)
            self.pad_start(n = padsize)
        self.increment((rg,rg),(np.sum(read.not_skipped_errors()), ~read.skips))

class QCovariate(Covariate):
    """
    A 2d covariate array in the read group and Q dimensions.
    """
    def __init__(self):
        super().__init__(shape = (0,0))

    def consume_read(self, read):
        rg = read.get_rg_int()

        #check boundaries
        largest_rg_idx = self.shape()[0] - 1
        if rg > 

class PosCovariate(Covariate):
    """
    A 3d covariate array in the read group, Q, and position dimensions.
    """
    def __init__(self):
        super().__init__(shape = (0,0,0))

class DinucCovariate(Covariate):
    """
    A 3d covariate array in the read group, Q, and dinucleotide dimensions.
    """
    def __init__(self):
        super().__init__(shape = (0,0,0))

class CovariateData():
    """
    A class for holding all the covariate data needed
    to recalibrate reads.

    Read-only Attributes

        Don't modify these or things will break

        * :attr:`_seqlen` - Sequence length
        * :attr:`_nrgs` - Number of read groups
        * :attr:`_maxscore` - Highest quality score

    Attributes

        * :attr:`rgcov` - RGCovariate
        * :attr:`qcov` - QCovariate
        * :attr:`poscov` - PosCovariate
        * :attr:`dinuccov` - DinucCovariate


    Methods

        * :meth:`consume_read` - Add covariates from the read to the proper data arrays

    """
    #nrgs will be available from the read class
    seqlen = 0

