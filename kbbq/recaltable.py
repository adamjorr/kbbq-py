#!usr/bin/env python3
import numpy as np
import pandas as pd

class GATKReport:
    """
    A base class representing a GATK Report file that contains many tables.

    Attributes

        * :attr:`header` - a header line (such as ``#:GATKReport.v1.1:5``)
        * :attr:`tables` - a list of GATKTable objects

    Methods

        * :meth:`write` - write the report to a file
        * :meth:`fromfile` - read the report from a file
    """

    def __init__(self, header, tables):
        """
        Initialize a GATKReport.

        :param str header: The header string
        :param list(str) tables: a list of tables.
        :ret: a representation of a report
        :rtype: :class:`.GATKReport`
        """
        self.header = header
        """A header line (such as ``#:GATKReport.v1.1:5``)"""
        self.tables = tables
        """A list of :class:`.GATKTable` objects."""

    @classmethod
    def fromfile(cls, filename):
        """
        Initialize the GATKReport from a file.

        :param str filename: The file to be read.
        :ret: a representation of a report
        :rtype: :class:`.GATKReport`
        """
        with open(filename) as fh:
            header = fh.readline()
            table_strings = fh.read().split('\n\n')
            tables = [GATKTable.fromstring(s) for s in table_strings if s != '']
            return cls(header, tables)
            

    def write(self, filename):
        """
        Write the report to filename.

        :param str filename: The file name to write.
        """
        with open(filename, 'w') as fh:
            fh.write(str(self))

    def __str__():
        """
        A string containing the report.

        Each table is separated by two newlines.
        """
        return str(self.header) + '\n' +
            '\n\n'.join([str(t) for t in self.tables] + [''])

class GATKTable:
    """
    A class representing a GATK report table.

    The intended method of interacting with the table is using the data
    attribute, which is a standard Pandas dataframe.

    Attributes

        * :attr:`format` - the format string representing the data in the table.
            This is like ``#:GATKTable:ncol:nrow:%f:%f:%f:;``.
            It is ``:`` delimited and ``;`` terminated.
        * :attr:`ncols` - The number of columns in the table
        * :attr:`nrows` - The number of rows in the table
        * :attr:`fmtstrings` - A list of the format strings for each column in the table
        * :attr:`name` - The entire line specifying the name of the table. It is
            formatted like ``#:GATKTable:title:subtitle``. ``subtitle`` can be
            an empty string.
        * :attr:`title` - The title of the table
        * :attr:`description` - The a description of the table. May be ``''``.
        * :attr:`data` - A Pandas dataframe containing the table data. Accessing this
            attribute is the primary way to interact with the data.

    Methods

        * :meth:`fromstring` - Initialize a :class:`GATKTable` from a string.
        * :meth:`parse_fmtstring` - Parse a fmtstring into a dictionary of types
        * :meth:`get_fmtstring` - Get a formatstring by inspecting the dataframe.
        * :meth:`get_titlestring` - Get the title string
        * :meth:`get_datastring` - Format the dataframe into a string
        * :meth:`get_nrows` - Get the number of rows in the dataframe.
        * :meth:`get_ncols` - Get the number of columns in the dataframe.
        * :meth:`write` - Write the table to a filehandle.
    """

    def __init__(self, title, description, data):
        """
        Initialize the table.

        TODO: refactor such that the doesn't store ncols, nrows, fmtstrings
        and uses a function that inspects the data to create the format line.
        Do the same with title and subtitle to recreate the name line.
        Header should be a part of data and not its own attribute.
        These ensure the attributes are all consistent even if the underlying
        data are manipulated by the user. (in short, move to getter functions)

        :param str title:
        :param str description:
        :param pandas.DataFrame data:

        """

        self.title = title
        """The title of the table."""
        self.description = description
        """A short description of the table."""
        self.data = data
        """
        A Pandas dataframe containing the table data. Accessing this
        attribute is the primary way to interact with the data.
        """
        self.typemap = {np.dtype(np.int) : '%d', np.dtype(np.float) : '%f', str : '%s', np.dtype(np.object) : '%s'}
        """
        A dictionary of `{type : str}` to determine the string for each type.

        Note that in the current implementation this only affects
        :meth:`get_fmtstring`.
        """

    @classmethod
    def fromstring(cls, tablestring):
        """
        Initialize the table from a string

        The first line of the table contains the a string describing the table
        format and is like ``#:GATKTable:ncol:nrow:%f:%f:%f:;``. It is ``:``
        delimited and ``;`` terminated. The last format string will be empty.
        The second line of the table specifies the name of the table. It is
        formatted like ``#:GATKTable:title:description``. :attr:`description`
        can be an empty string. It is ``:`` delimited. The third line contains
        the header. The header and row data are delimited by two or more spaces.
        The number of spaces are determined by padding each value in the row
        such that the the character width of the column is constant and can
        represent all the  data in the column. Strings (including column
        headings) are left justified and numeric types are right justified.
        """

        rows = tablestring.splitlines()
        title, description = rows[1].split(':')[2:4]
        header = rows[2].split()
        typedict = cls.parse_fmtstring(header, rows[0])
        strdata = [s.split() for s in rows[3:]]
        d = dict(zip(header, zip(*strdata))) #dictionary {colname : coldata}
        data = pd.DataFrame(d).astype(typedict)
        return cls(title, description, data)

    @staticmethod
    def parse_fmtstring(header, fmtstring):
        """
        Turn a fmtstring into a dictionary of `{col_title : type}`

        :param list(str) header:
        :param str fmtstring:
        :return: A dictionary mapping `{col_title : type}`
        :rtype: dict(:class:`str`, :obj:`type`)
        """
        splitfmt = fmtstring.split(':')
        fmtstrings = splitfmt[4:-1]
        typedict = {}
        for i,h in enumerate(header):
            f = fmtstrings[i]
            if f.endswith('d'):
                t = np.int64
            elif f.endswith('f'):
                t = np.float64
            elif f.endswith('s'):
                t = str #pandas will convert this to 'object'
            else:
                t = None
            if t is not None:
                typedict[h] = t
        return typedict

    def get_fmtstring(self):
        """
        Get the format string by inspecting the dataframe.

        Uses the dtypes, number of rows, and number of columns in the data
        frame along with :attr:`typemap` to create the format string, which is
        formatted like ``#:GATKTable:ncol:nrow:%f:%f:%f:;``.

        :return: The format string
        :rtype: str

        """
        fmtlist = ['#', 'GATKTable', str(self.get_ncols()), str(self.get_nrows())]
        types = self.data.dtypes
        for t in types:
            fmtlist.append(self.typemap[t])
        fmtlist.append(';')
        return ':'.join(fmtlist)

    def get_titlestring(self):
        """
        Get the titlestring.

        Uses the :attr:`title` and :attr:`description` to create the title
        line. It is formatted like ``#:GATKTable:title:subtitle.``

        :return: the title string
        :rtype: str
        """
        titlelist = ['#', 'GATKTable', self.title, self.description]
        return ':'.join(titlelist)

    def get_datastring(self):
        """
        Get the data string. Includes the header line and the data in
        the table.

        :return: the data table as a formatted string
        :rtype: str
        """
        datawidths = self.data.apply(lambda x: x.str.len().max()).to_numpy()
        headwidths = self.data.columns.str.len().to_numpy()
        header = self.data.columns.to_list()
        colwidths = np.maximum(datawidths, headwidths)
        fmstrings = [self.typemap.get(t) for t in self.data.dtypes]
        datastr = '  '.join([h.ljust(colwidths[i]) for i,h in enumerate(header)])
        for row in self.data.itertuples(index = False):
            formatted = [getattr(row,header[i]).ljust(colwidths[i]) if \
                f == '%s' else \
                (f % float(getattr(row,header[i]))).rjust(colwidths[i]) \
                for i,f in enumerate(fmtstrings)]
            datastr = datastr + '\n' + '  '.join(formatted)
        return datastr

    def get_nrows(self):
        """
        Get the number of rows in the dataframe.

        :return: The number of rows in the dataframe.
        :rtype: int
        """
        return self.data.shape[0]

    def get_ncols(self):
        """
        Get the number of columns in the dataframe.

        :return: The number of columns in the dataframe.
        :rtype: int
        """
        return self.data.shape[1]

    def write(self, filehandle):
        """
        Write the table to a filehandle.

        Adds a newline to the end that isn't present in :meth:`__str__`
        :param filehandle: The file object to write the table to.
        :type filehandle: File Object
        :return: The file object's return (usually the number of characters written)
        :rtype: Usually int

        """
        return filehandle.write(str(self) + '\n')

    def __str__(self):
        return self.get_fmtstring() + '\n' + self.get_titlestring() + '\n' + get_datastring()

    def __repr__(self):
        return self.get_fmtstring() + '\n' + self.get_titlestring() + '\n' + repr(self.data)

class RecalibrationReport(GATKReport):
    """
    A class representing a GATK Recalibration Report used for Base Quality 
    Score Recalibration.

    It extends the :class:`GATKReport`, so the primary
    method of interacting with the table is the :attr:`tables` attribute.
    The primary difference is the :func:`__init__` function will set indices
    and data types as they should be for a recalibration report.

    * Table :attr:`tables[0]` contains the arguments used to call the program
        that generated the table.
    * Table :attr:`tables[1]` contains the quanization map
    * Table :attr:`tables[2]` contains read group error information.
    * Table :attr:`tables[3]` contains error information subset by read group,
        then reported quality score.
    * Table :attr:`tables[4]` contains error information subset by read group,
        then reported quality score, then covariate name and value. Both
        Context and Cycle are valid for the covariate name.

    See the table below for a schema.

    =====  =================  ===========================================
    Table  Number of Columns  Variables
    =====  =================  ===========================================
    0      2                  Argument, Value
    1      3                  QualityScore, Count, QuantizedScore
    2      6                  ReadGroup, QualityScore, EventType,
                              EmpiricalQuality, Observations, Errors
    3      6                  ReadGroup, QualityScore, EventType,
                              EmpiricalQuality, Observations, Errors
    4      8                  ReadGroup, QualityScore, CovariateValue,
                              CovariateName, EventType, EmpiricalQuality,
                              Observations, Errors
    =====  =================  ===========================================

    For current GATK parameters, EventType will be ``M`` for mismatch though
    ``I`` for insertion and ``D`` for deletion are possible. Valid values for
    CovariateName are ``Context`` and ``Cycle``.

    """

    def __init__(self, filename):
        """
        Initialize the recalibration report from a file.

        This initializes the tables with some data wrangling to set indices
        and data types properly. If you don't want this, use a 
        :class:`.GATKReport` instead.
        """
        super().__init__(filename)
        #self.data[0] is argument / value
        #self.data[1] is quantization map
        #self.data[2] is RG
        #self.data[3] is RG / reportedqual
        #self.data[4] is RG / reportedqual / covariate (Cycle OR Context)
        self.tables[0].data = self.tables[0].data.set_index('Argument')
        self.tables[1].data = self.tables[1].data.astype({'QualityScore' : np.int, 'Count' : np.longlong, 'QuantizedScore' : np.int })
        self.tables[1].data = self.tables[1].data.set_index('QualityScore')
        #self.tables[2].data = self.tables[2].data.astype(typer, errors = 'ignore')
        self.tables[2].data = self.tables[2].data.set_index('ReadGroup')
        #typer.pop('EstimatedQReported')
        #typer.update({'QualityScore' : np.int_})

        typer = {'QualityScore' : np.int_}
        self.tables[3].data = self.tables[3].data.astype(typer)
        self.tables[3].data = self.tables[3].data.set_index(['ReadGroup','QualityScore'])
        self.tables[4].data = self.tables[4].data.astype(typer)
        self.tables[4].data = self.tables[4].data.set_index(['ReadGroup','QualityScore','CovariateName','CovariateValue'])


