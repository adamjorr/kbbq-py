"""
Class and functions for read data.

This module includes the :class:`ReadData` class and a few utility functions for 
working with BAM reads of type :class:`pysam.AlignedSegment`.

Classes

    * :class:`ReadData`

Functions

    * :func:`bamread_get_oq`
    * :func:`bamread_get_quals`

"""

class ReadData():
    """
    A class that represents some minimal information from a sequencing read.

    Since this won't hold as much data as a BAM read, if you want
    to recalibrate a read you should turn it into a ReadData object,
    then update the original read with the new information; otherwise
    you'll lose annotations. This class also manages read group information.
    You should instantiate the object from one of the class methods instead
    of directly instantiating it if possible.

    Class Attributes

        * :attr:`rg_to_pu` - Dict of read group id's -> platform unit
        * :attr:`rg_to_int` - Dict of read group id's -> int
        * :attr:`numrgs` - int number of read groups

    Instance Attributes

        * :attr:`seq` - Numpy array of characters representing the sequence
        * :attr:`qual` - Numpy array of int representing the quality score
        * :attr:`skips` - Numpy array of bools representing sites to skip
        * :attr:`name` - string representing the name of the read
        * :attr:`rg` - int representing the read group the read belongs to
        * :attr:`second` - bool representing whether the read is 2nd in pair
        * :attr:`errors` - Numpy array of bools representing whether the base is an error

    Class Methods

        * :meth:`from_bamread` - instantiate a ReadData object from a BAM read 
        * :meth:`from_fastq` -  instantiate a ReadData object from a fastq read
        * :meth:`load_rgs_from_bamfile` - load read groups into the class from a bam file object

    Instance Methods

        * :meth:`str_qual` - Get the quality score as a list of chars
        * :meth:`canonical_name` - Return the name with a '/1' or '/2' suffix
        * :meth:`get_rg_int` - Get the read group index using the rg_to_int dictionary
        * :meth:`get_pu` - Get the PU from the read group and rg_to_pu dictionary.
        * :meth:`not_skipped_errors` - Return a logical and of ~:attr:`skips` and :attr:`errors`

    Note that from https://docs.python.org/3/library/stdtypes.html#dict , iteration
    order is guaranteed to be in insertion order. Thus we are OK saving the rg as an
    int on the fly so long as we don't remove any read groups from the rg_to_pu dictionary.
    So don't do that! In fact, we should probably remove __del__() and pop() from the dicts...
    """

    rg_to_pu = dict()
    rg_to_int = dict()
    numrgs = 0

    def __init__(self, seq, qual, skips, name, rg, second, errors):
        self.seq = seq
        self.qual = qual
        self.skips = skips
        self.name = name
        self.rg = rg
        if rg not in rg_to_pu:
            #if it hasn't been preloaded,
            #we create a new PU identical to the rg
            #and load it
            rg_to_pu[rg] = rg
            rg_to_int[rg] = numrgs
            numrgs = numrgs + 1
        self.second = second
        self.errors = errors

    @classmethod
    def from_bamread(cls, bamread, use_oq = False):
        """
        ReadData factory that instantiates a ReadData object from a pysam AlignedSegment.

        Skips and errors are initialized to be empty bool arrays. If you need to use
        these attributes, make sure you set them the way you need to for your use case.

        If :code:`use_oq` is :code:`True`, use the OQ tag for base quality. 
        The default is to use regular quality scores.

        If the read group hasn't been loaded in the dictionary, it will be registered
        with a PU equal to the value of :code:`rg`. To load the dictionary, call
        :meth:`load_rgs_from_bamfile` first.

        :param bamread: read to get data from
        :type bamread: :class:`pysam.AlignedSegment`
        :param bool use_oq: use the OQ tag for quality scores
        :return: a read object
        :rtype: :class:`ReadData`
        """
        seq = np.array(list(bamread.query_sequence), dtype = np.unicode)
        seqlen = len(seq)
        return cls(
            seq = seq,
            qual = bamread_get_quals(bamread, use_oq),
            skips = np.zeros(seqlen, dtype = np.bool),
            name = bamread.query_name,
            rg = bamread.get_tag('RG'),
            second = bamread.is_read2,
            errors = np.zeros(seqlen, dtype = np.bool),
            )

    @classmethod
    def from_fastq(cls, fastqread, rg = None, second = None, namedelimiter = '_'):
        """
        ReadData factory that instantiates a ReadData object from a pysam FastqProxy.

        Skips and errors are initialized to be empty bool arrays. If you need to use them,
        make sure you set them the way you need to for your use case. The read name will
        be set to the first field of the delimited read name, minus any trailing :code:`/1` or :code:`/2`
        if they exist.

        If :code:`rg` is None (the default), we will attempt to infer the read group id from the
        read name. To infer rg, the read group must be
        in its own field with the field beginning with :code:`RG:`, such as :code:`RG:example`,
        with fields delimited by :code:`namedelimiter`. When multiple fields are found that
        begin with :code:`RG:`, the last is chosen to be the true ID.

        If the read group (either inferred or provided explicitly) hasn't been
        loaded in the dictionary, it will be registered with a PU equal to the value of :code:`rg`.

        If second is None (the default) we will attempt to infer if the read
        is second in pair based on the name.
        To infer second, the first field of the read name must end with :code:`/2`.
        If the last 2 characters of the first field are not :code:`/2`, second will be inferred to be false.

        :param fastqread: The fastq read to get data from
        :type fastqread: :class:`pysam.FastqProxy`
        :param str rg: the read group the read belongs to
        :param bool second: whether the read is 2nd in pair
        :param str namedelimiter: the delimiter for parsing the read name
        :return: a read object
        :rtype: :class:`ReadData`
        """
        seq = np.array(list(read.sequence), dtype = np.unicode)
        seqlen = len(seq)
        splitname = read.name.split(sep = namedelimiter)
        if rg is None:
            rg = [f[3:] for f in splitname if f[0:3] == 'RG:'][-1]
        if second is None:
            second = (splitname[0][-2:] == '/2')
        if splitname[0].endswith(('/1','/2')):
            splitname[0] = splitname[0][:-2]

        return cls(
            seq = seq,
            qual = np.array(read.get_quality_array(), dtype = np.int),
            skips = np.zeros(seqlen, dtype = np.bool),
            name = splitname[0],
            rg = rg,
            second = second,
            errors = np.zeros(seqlen, dtype = np.bool),
            )

    @classmethod
    def load_rgs_from_bamfile(cls, bamfileobj):
        """
        Load read group IDs and PUs from the header of the
        pysam bamfileobj.

        Recalibration is done on a PU basis, so when building
        a model using reads in a BAM, the PU for each read is
        obtained by looking up the PU associated with its read
        group. The dictionary that controls these lookups is a
        class variable. To store the RG data in an int form, an
        rg-to-int dict is also loaded. This should be equal to
        :code:`dict(zip(rg_to_pu, range(len(rg_to_pu))))`

        :param bamfileobj: The opened bam file
        :type bamfileobj: :class:`pysam.AlignmentFile`
        """
        for rg in bamfileobj.header.as_dict()['RG']:
            rg_to_pu[rg['ID']] = rg['PU']
            rg_to_int[rg['ID']] = numrgs
            numrgs = numrgs + 1

    def str_qual(self, offset = 33):
        """
        Get the quality of this read as a list of characters.

        offset (default 33) will be added to the ASCII value
        of each character.

        :param int offset: Offset to add
        :return: read quality
        :rtype: list(chr)
        """
        return list((self.qual + offset).astype(np.uint32).view('U1'))

    def canonical_name(self):
        """
        The name with an added suffix based on whether the read
        is firstinpair or not.

        If the read has its second flag set to True, :code:`/2` is added
        to the end. Otherwise, :code:`/1` is added.

        :return: read name with a suffix
        :rtype: str
        """
        suffix = ("/2" if self.second else "/1")
        return self.name + suffix

    def get_rg_int(self):
        """
        Return the RG as an int suitable for indexing rather
        than as the actual string RG ID.

        :return: read group index
        :rtype: int
        """
        return rg_to_int[self.rg]

    def get_pu(self):
        """
        Return the PU from the rg_to_int dict.

        :return: The Read Group's PU tag
        :rtype: str
        """
        return rg_to_int[self.rg]

    def not_skipped_errors(self):
        """
        Return a logical and of ~:attr:`skips` and :attr:`errors`

        :return: array of valid errors
        :rtype: :class:`np.ndarray`(bool)
        """
        return np.logical_and(self.errors, ~self.skips)

    def __len__(self):
        """
        Return sequence length.

        :return: sequence length
        :rtype: int
        """
        return len(self.seq)

def bamread_get_oq(read, offset = 33):
    """
    Get the OQ of the given bamread as an array of int.

    offset (default 33) will be subtracted from the ASCII value
    of each character.

    :param read: Read to get quals from
    :type read: :class:`pysam.AlignedSegment`
    :param int offset: Offset to subtract
    :return: quality scores
    :rtype: :class:`numpy.ndarray` of int

    """
    oq = np.array(list(read.get_tag('OQ')), dtype = np.unicode)
    quals = np.array(oq.view(np.uint32) - offset, dtype = np.uint32)
    return quals

def bamread_get_quals(read, use_oq = False):
    """
    Return the qualities of a pysam bam read as an array.

    If use_oq = True, use the OQ tag for base quality.
    The default is to use regular quality scores.
    
    :param read: Read to get quals from
    :type read: :class:`pysam.AlignedSegment`
    :param bool use_oq: Use OQ tag for quality scores
    :return: quality scores
    :rtype: :py:class:`numpy.ndarray` of int

    """
    if use_oq:
        return bamread_get_oq(read)
    else:
        return np.array(read.query_qualities, dtype = np.int)
