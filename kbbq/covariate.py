"""
A data class and helper functions for read covariates.
"""
from . import read
from . import compare_reads

def pad_axis(array, axis, n):
    """
    Pad the axis of :class:`np.ndarray` with n zeros.
    """
    return np.append(array, np.zeros(array.shape[0:axis] + (n,) + array.shape[axis+1:]), axis = axis)

class Covariate():
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
        
    def pad_axis(self, axis, n = 1):
        """
        Expand axis dimension by n.
        """
        self.errors = pad_axis(self.errors, axis = axis, n = n)
        self.total = pad_axis(self.total, axis = axis, n = n)

    def pad_axis_to_fit(self, axis, idx):
        """
        Expand the axis until idx is a valid idx in that dimension.

        Does nothing if idx is already valid. idx should be a single
        scalar value, not a slice.
        """
        axislen = self.shape()[axis]
        if idx < -axislen or idx >= axislen:
            #out of bounds
            if idx < 0:
                padsize = -idx - axislen 
            else:
                padsize = idx - axislen + 1
            self.pad_axis(axis = axis, n = padsize)

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
        rge, rgv = read.get_rg_errors()
        self.pad_axis_to_fit(axis = 0, idx = rgv[0])
        self.increment((rge,rgv))
        return rge, rgv

    def num_rgs(self):
        return self.shape()[0]

class QCovariate(Covariate):
    """
    A 2d covariate array in the read group and Q dimensions.
    """
    def __init__(self):
        self.rgcov = RGCovariate()
        super().__init__(shape = (0,0))

    def consume_read(self, read):
        rge, rgv = self.rgcov.consume_read(read)
        self.pad_axis_to_fit(axis = 0, idx = self.rgcov.num_rgs() - 1)

        qe, qv = read.get_q_errors(read)
        try:
            self.increment(idx = ((rge,qe),(rge,qv)))
        except IndexError:
            #this way we only need to find the max if assignment fails
            self.pad_axis_to_fit(axis = 1, idx = np.amax(qv))
            self.increment(idx = ((rge,qe),(rge,qv)))
        return (rge, rgv), (qe, qv)

    def num_qs(self):
        return self.shape()[1]

class CycleCovariate(Covariate):
    """
    A 3d covariate array in the read group, Q, and cycle dimensions.
    """
    def __init__(self):
        super().__init__(shape = (0,0,0))

class DinucCovariate(Covariate):
    """
    A 3d covariate array in the read group, Q, and dinucleotide dimensions.
    """
    def __init__(self):
        super().__init__(shape = (0,0,len(compare_reads.Dinucleotide.dinucs)))

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
        * :meth:`get_num_rgs` - Return number of rgs in the rg covariate
        * :meth:`get_num_qs` - Return number of qs in the q covariate
        * :meth:`get_num_cycles` - Return number of cycles in the cycle covariate
        * :meth:`get_num_dinucs` - Return number of dinucs in the dinuc covariate
    """
    def __init__(self):
        self.qcov = QCovariate()
        self.cyclecov = CycleCovariate()
        self.dinuccov = DinucCovariate()

    def consume_read(self, read):
        (rge, rgv), (qe, qv) = self.qcov.consume_read(read)
        ce, cv = read.get_cycle_errors()
        de, dv = read.get_dinuc_errors()

        num_rgs = self.qcov.rgcov.num_rgs()
        num_qs = self.qcov.num_qs()
        for cov in self.cyclecov, self.dinuccov:
            cov.pad_axis_to_fit(axis = 0, idx = num_rgs - 1)
            cov.pad_axis_to_fit(axis = 1, idx = num_qs - 1)

        self.cyclecov.pad_axis_to_fit(axis = 2, idx = 2 * len(read) - 1)
        self.cyclecov.increment(idx = ((rge,qe,ce),(rgv,qv,cv)))
        self.dinuccov.increment(idx = ((rge,qe,de),(rgv,qv,dv)))
