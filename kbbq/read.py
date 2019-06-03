"""
Class and functions for read data.
"""

class ReadData():
    """
    A class that represents some minimal informatino from a sequencing read.

    Since this won't hold as much data as a BAM read, if you want
    to recalibrate a read you should turn it into a ReadData object,
    then update the original read with the new information; otherwise
    you'll lose annotations. This class also manages read group information.
    You should instantiate the object from one of the class methods instead
    of directly instantiating it if possible.

    Class Attributes

        * attr:`rgs` - List of observed read groups
        * attr:`pus` - List of platform unit tags associated with the read groups

    Instance Attributes

        * attr:`seq` - Numpy array of characters representing the sequence
        * attr:`qual` - Numpy array of int representing the quality score
        * attr:`skips` - Numpy array of bools representing sites to skip
        * attr:`name` - string representing the name of the read
        * attr:`rg` - int representing the read group the read belongs to
        * attr:`second` - bool representing whether the read is 2nd in pair

    Class Methods

        * meth:`from_bamread` - instantiate a ReadData object from a BAM read 
        * meth:`from_fastq` -  instantiate a ReadData object from a fastq read
        * meth:`get_rg_to_int` - get a dict mapping read group ID to an int
        * meth:`get_rg_to_pu` - get a dict mapping read group ID to a PU tag

    Instance Methods

        * meth:`str_qual` - Get the quality score as a list of chars

    """

    rgs = []
    pus = []

    def __init__(self, seq, qual, skips, name, rg, pu, second):
        self.seq = seq
        self.qual = qual
        self.skips = skips
        self.name = name
        self.rg = rg
        self.pu = pu
        if rg not in rgs:
            #this may be slow...
            rgs.append(rg)
            pus.append(pu)
        self.second = second

    @classmethod
    def from_bamread(cls, bamread):
        """
        TODO
        """
        pass

    @classmethod
    def from_fastq(cls, fastqread):
        """
        TODO
        """
        pass

    @classmethod
    def get_rg_to_int(cls):
        return {r:i for i,r in enumerate(rgs)}

    @classmethod
    def get_rg_to_pu(cls):
        return dict(zip(rgs, pus))

    def str_qual(self, offset = 33):
        """
        Get the quality of this read as a list of characters.
        """
        return list((self.qual + offset).astype(np.uint32).view('U1'))


def bamread_get_oq(read, offset = 33):
    """
    Get the OQ of the given bamread.
    """
    oq = np.array(list(read.get_tag('OQ')), dtype = np.unicode)
    quals = np.array(oq.view(np.uint32) - offset, dtype = np.uint32)
    return quals

def bamread_get_quals(read, use_oq = False):
    """
    Return the qualities of the read as an np.array.

    Use the OQ tag for read quals if use_oq = True.
    """
    if use_oq:
        return bamread_get_oq(read)
    else:
        return np.array(read.query_qualities, dtype = np.int)
