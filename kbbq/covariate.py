"""
A data class and helper functions for read covariates.

This module includes the :class:`Covariate` class and subclasses,
as well as the :func:`pad_axis` helper function for appending
to the end of an arbitrary axis of an :class:`numpy.ndarray`.

Classes
    * :class:`Covariate`
    * :class:`RGCovariate`
    * :class:`QCovariate`
    * :class:`CycleCovariate`
    * :class:`DinucCovariate`

Functions
    * :func:`pad_axis`
"""
from . import read
from . import compare_reads

def pad_axis(array, axis, n):
    """
    Pad the axis of :class:`np.ndarray` with n zeros.

    :param int axis: axis to pad
    :param int n: size of pad to apply
    :return: padded array
    :rtype: :class:`numpy.npdarray`
    """
    return np.append(array, np.zeros(array.shape[0:axis] + (n,) + array.shape[axis+1:]), axis = axis)

class Covariate():
    """
    A base class for covariates.

    Holds 2 :class:`np.ndarray`, :attr:`errors` for number of errors and
    :attr:`total` for the total number of observations of this covariate.
    The indices of each array are they keys; ie. the number of times
    covariate 2 was observed can be found with :code:`self.total[2]`.

    .. warning::
       Don't assign directly to :attr:`errors` or :attr:`total`.
       :meth:`shape` assumes both arrays have the same shape and
       only returns the shape of :attr:`total`. If the arrays'
       shapes become different somehow, :meth:`__getitem__`
       and :meth:`__setitem__` may start throwing errors unexpectedly.
       If you want to use two arrays that aren't guaranteed to have
       the same size, use two arrays. This isn't the right class
       for that use case.

    :meth:`__getitem__` will pass the key to both arrays and return a tuple
    of (errors, total). As a pneumonic, you can remember this as a fraction
    from left to right, as total should always be larger than errors.

    :meth:`__setitem__` will also pass the key to both arrays, and will expect
    value to be a tuple of values, setting the first value to errors and the
    second to total.

    Attributes

        * :attr:`errors` - :class:`numpy.ndarray` of the number of errors observed
        * :attr:`total` - :class:`numpy.ndarray` of the total number of observations


    Methods

        * :meth:`__init__` - Initialize arrays with desired shape.
        * :meth:`pad_axis` - Pad specified axis of both arrays
        * :meth:`pad_axis_to_fit` - Pad specified axis so the given index is valid
        * :meth:`increment` - Increment the given indices of each array.
        * :meth:`shape` - return the shape of the arrays.
        * :meth:`__getitem__` - get number of errors and total observations at the given index
        * :meth:`__setitem__` - set errors and observations at the given index

    """
    def __init__(self, shape = 0):
        """
        Initialize the arrays with the given shapes.

        :param shape tuple: desired shape of error and total arrays
        :return: empty object
        :rtype: :class:`Covariate`
        """
        self.errors = np.zeros(shape, dtype = np.int)
        """:class:`np.ndarray` (int) holding number of errors observed for each key."""
        self.total = np.zeros(shape, dtype = np.int)
        """:class:`np.ndarray` (int) holding number of observations for each key."""
        
    def pad_axis(self, axis, n = 1):
        """
        Expand axis dimension by n.

        :param int axis: axis to pad
        :param int n: size of pad to apply
        """
        self.errors = pad_axis(self.errors, axis = axis, n = n)
        self.total = pad_axis(self.total, axis = axis, n = n)

    def pad_axis_to_fit(self, axis, idx):
        """
        Expand the axis until idx is a valid idx in that dimension.

        Does nothing if idx is already valid. idx should be a single
        scalar value, not a slice.

        :param int axis: axis to pad
        :param int idx: index to ensure is valid
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
        Increment the given indices by value.

        :code:`idx` should be a tuple of values interpretable by :func:`np.add`.

        The first indices given will increment the error array.
        The second indices given will increment the total array.

        Note that this is implemented with :func:`np.add.at`.
        For example, if :code:`idx[0]` is :code:`[2,2]`, :attr:`errors` [2]
        will be incremented twice. This is equivalent to

        .. highlight:: python

        .. code-block:: python

           errors[2] += 1
           errors[2] += 1
        
        If :code:`value[0]` is 2, :attr:`errors` will be instead
        incremented by 4 and is equivalent to

        .. code-block:: python

           errors[2] += 2
           errors[2] += 2

        :param idx: indices to increment array at
        :type idx: tuple( index , index )
        :param value: value to add for each index given.
        :type value: typle( int, int )
        """
        np.add.at(self.errors, idx[0], value[0])
        np.add.at(self.total, idx[1], value[1])

    def shape(self):
        """
        Return the shape of the underlying arrays.

        :return: array shape
        :rtype: tuple(int)
        """
        return self.total.shape

    def __getitem__(self, key):
        """
        Pass key to each :class:`np.ndarray` and return the result as a tuple.

        The result will be in order of (errors, total). As a pneumonic,
        you can remember this as a fraction from left to right,
        as total should always be larger than errors.

        :param slice key: indices to get
        :return: number of errors and number of total observations
        :rtype: tuple(int, int) or tuple( :class:`numpy.ndarray` , :class:`numpy.ndarray` )
        """
        return (self.errors[key], self.total[key])

    def __setitem__(self, key, value):
        """
        Pass key to each :class:`np.ndarray` and assign each value.

        Value must be a tuple of values to assign. The first value
        is assigned to errors and the second to total.

        :param slice key: indices to set
        :param tuple(int,int) value: values to set
        """
        self.errors[key] = value[0]
        self.total[key] = value[1]

class RGCovariate(Covariate):
    """
    A 1d covariate array in the read group dimension.

    Methods

        * :meth:`__init__` - Initialize arrays
        * :meth:`consume_read` - increment the covariates given a :class:`kbbq.read.ReadData`
        * :meth:`num_rgs` - return number of rgs observed so far

    """
    def __init__(self):
        """
        Initialize a 1D read group array.

        :return: empty object
        :rtype: :class:`Covariate`
        """
        super().__init__(shape = 0)

    def consume_read(self, read):
        """
        Increment the rg covariates given a :class:`kbbq.read.ReadData`

        :param read: read to take covariate values from
        :type read: :class:`kbbq.read.ReadData`
        :return: erroneous and valid read group values
        :rtype: tuple( :class:`numpy.ndarray` , :class:`numpy.ndarray` )
        """
        rge, rgv = read.get_rg_errors()
        self.pad_axis_to_fit(axis = 0, idx = rgv[0])
        self.increment((rge,rgv))
        return rge, rgv

    def num_rgs(self):
        """
        Return size needed to hold largest rg observed so far.

        :return: largest rg + 1
        :type: int
        """
        return self.shape()[0]

class QCovariate(Covariate):
    """
    A 2d covariate array in the read group and Q dimensions.

    This class holds a :class:`RGCovariate` to efficiently increment
    both objects while iterating through reads.

    Attributes

        * :attr:`rgcov` - RG covariates

    Methods

        * :meth:`__init__` - Initialize arrays
        * :meth:`consume_read` - increment the rg and q covariates given a :class:`kbbq.read.ReadData`
        * :meth:`num_qs` - return number of qs observed so far

    """
    def __init__(self):
        """
        Initialize an empty :class:`RGCovariate` and initialize the empty q arrays.
        """
        self.rgcov = RGCovariate()
        """RG covariates"""
        super().__init__(shape = (0,0))

    def consume_read(self, read):
        """
        Increment the rg and q covariates given a :class:`kbbq.read.ReadData`

        :param read: read to take covariate values from
        :type read: :class:`kbbq.read.ReadData`
        :return: erroneous and valid read group values and erroneous and valid q values
        :rtype: (( :class:`numpy.ndarray` , :class:`numpy.ndarray` ),( :class:`numpy.ndarray` , :class:`numpy.ndarray` )) 
        """
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
        """
        Return size needed to hold largest q observed so far.

        :return: largest q + 1
        :type: int
        """
        return self.shape()[1]

class CycleCovariate(Covariate):
    """
    A 3d covariate array in the read group, Q, and cycle dimensions.

    :meth:`pad_axis` is overriden here since the way cycles are stored
    means appending to the end of the array changes the location of the
    data.

    Methods:

        * :meth:`__init__` - Initialize empty arrays
        * :meth:`pad_axis` - Pad given axis.

    """
    def __init__(self):
        """
        Initialize empty arrays.
        """
        super().__init__(shape = (0,0,0))

    def pad_axis(self, axis, n = 1):
        """
        Expand the axis dimension by n.

        This is overriden for the CycleCovariate since we use one array
        to handle both positive and negative cycle values. Appending
        0's to the end would ruin the array since negative indices are
        significant. Thus, we have to wrangle the existing data a bit.

        :param int axis: axis to pad
        :param int n: size of pad to apply
        """
        if not (axis == 2 or axis == -1):
            super().pad_axis(self, axis = axis, n = n)
        else:
            if n % 2 != 0:
                raise ValueError('n should be even for the 2nd axis ' +
                    'of a CycleCovariate. n = {} was given.'.format(n))
            oldlen = self.shape()[2]
            newlen = oldlen + n
            newerrors = np.zeros(self.shape()[0:2] + (newlen,), dtype = np.int)
            newtotal = np.zeros(self.shape()[0:2] + (newlen,), dtype = np.int)
            newerrors[...,0:oldlen/2], newtotal[...,0:oldlen/2] = self[...,0:oldlen/2]
            newerrors[...,-oldlen/2:], newtotal[...,-oldlen/2:] = self[...,-oldlen/2:]
            self.errors = newerrors
            self.total = newtotal

    def num_cycles(self):
        """
        Return number of cycles in the cycle covariate

        This is the half the length of the cycle vector;
        there are twice as many cycles since positive and negative
        cycles are treated differently.

        :return: half the length of cycle vector
        :rtype: int
        """
        return self.shape()[-1] / 2

class DinucCovariate(Covariate):
    """
    A 3d covariate array in the read group, Q, and dinucleotide dimensions.
    """
    def __init__(self):
        """
        Initialize empty arrays.
        """
        super().__init__(shape = (0,0,len(compare_reads.Dinucleotide.dinucs)))

    def num_dinucs(self):
        """
        Return length of dinucleotide array.

        :return: number of dinucleotides
        :rtype: int
        """
        return self.shape()[-1]

class CovariateData():
    """
    A class for holding all the covariate data needed
    to recalibrate reads.

    Attributes

        * :attr:`rgcov` - RGCovariate
        * :attr:`qcov` - QCovariate
        * :attr:`poscov` - PosCovariate
        * :attr:`dinuccov` - DinucCovariate


    Methods - TODO

        * :meth:`consume_read` - Add covariates from the read to the proper data arrays
        * :meth:`get_num_rgs` - Return number of rgs in the rg covariate
        * :meth:`get_num_qs` - Return number of qs in the q covariate
        * :meth:`get_num_cycles` - Return number of cycles in the cycle covariate
        * :meth:`get_num_dinucs` - Return number of dinucs in the dinuc covariate
    """
    def __init__(self):
        """
        Initialize covariate arrays.
        """
        self.qcov = QCovariate()
        """Q covariates"""
        self.cyclecov = CycleCovariate()
        """Cycle covariates"""
        self.dinuccov = DinucCovariate()
        """Dinucleotide covariates"""

    def consume_read(self, read):
        """
        Add covariates from the read to the proper data arrays.

        :param read: the read to take data from
        :type read: :class:`kbbq.read.ReadData`
        """
        (rge, rgv), (qe, qv) = self.qcov.consume_read(read)
        ce, cv = read.get_cycle_errors()
        de, dv = read.get_dinuc_errors()

        num_rgs = self.get_num_rgs()
        num_qs = self.get_num_qs()
        for cov in self.cyclecov, self.dinuccov:
            cov.pad_axis_to_fit(axis = 0, idx = num_rgs - 1)
            cov.pad_axis_to_fit(axis = 1, idx = num_qs - 1)

        self.cyclecov.pad_axis_to_fit(axis = 2, idx = 2 * len(read) - 1)
        self.cyclecov.increment(idx = ((rge,qe,ce),(rgv,qv,cv)))
        self.dinuccov.increment(idx = ((rge,qe,de),(rgv,qv,dv)))

    def get_num_rgs(self):
        """
        Return number of rgs in the rg covariate

        :return: length of rg vector
        :rtype: int
        """
        return self.qcov.rgcov.num_rgs()

    def get_num_qs(self):
        """
        Return number of qs in the q covariate

        :return: length of q vector
        :rtype: int
        """
        return self.qcov.num_qs()

    def get_num_cycles(self):
        """
        Return number of cycles in the cycle covariate

        This is the half the length of the cycle vector;
        there are twice as many cycles since positive and negative
        cycles are treated differently.

        :return: half the length of cycle vector
        :rtype: int
        """
        return self.cyclecov.num_cycles()

    def get_num_dinucs(self):
        """
        Return number of dinucs in the dinuc covariate.

        :return: length of the dinuc vector
        :rtype: int
        """
        return self.dinuccov.num_dinucs()

