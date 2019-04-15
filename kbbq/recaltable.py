#!usr/bin/env python3
import numpy as np
import pandas as pd

class GATKReport:
    """
    A base class representing a GATK Report file that contains many tables.

    Attributes

        * :attr:`tables` - a list of GATKTable objects
        * :attr:`version` - the GATKReport format version. The default is
            ``'1.1'``.

    Methods

        * :meth:`write` - write the report to a file
        * :meth:`fromfile` - read the report from a file
        * :meth:`get_headerstring` - get the header line.
            Does not include a newline.
    """

    def __init__(self, tables, version = '1.1'):
        """
        Initialize a GATKReport.

        This module was written by reverse engineering a version 1.1 report.
        If you have a table with a newer version, please file a bug report,
        especially if this module fails to parse it.

        :param list(GATKTable) tables: a list of tables.
        :param str version: The report version as a string.
        :returns: a representation of a report
        :rtype: :class:`.GATKReport`

        """
        self.tables = tables
        """A list of :class:`.GATKTable` objects."""
        self.version = version
        """The report version identifier (such as ``'1.1'``)"""

    @classmethod
    def fromfile(cls, filename):
        """
        Initialize the GATKReport from a file.

        :param str filename: The file to be read.
        :returns: a representation of a report
        :rtype: :class:`.GATKReport`
        """
        with open(filename) as fh:
            fullheader = fh.readline()
            _, version, ntables = fullheader.strip().split(':')
            version = version.split(sep = 'v', maxsplit = 1)[-1]
            table_strings = fh.read().split('\n\n')
            tables = [GATKTable.fromstring(s) for s in table_strings if s != '']
            if (len(tables) != int(ntables)):
                raise ValueError(f"Malformed or truncated file {filename}.\
                    The header ({fullheader}) implies there should be\
                    {ntables} tables in this report, but we only found\
                    {len(tables)}. If you are sure the file is intact,\
                    modify the header line and try again.")
        return cls(tables, version)

    def get_headerstring(self):
        """
        Get the header line.

        The header line identifies the file as a report, shows the version,
        and shows the number of tables that should be present
        (such as ``#:GATKReport.v1.1:5``). This method returns a string
        that does not include a trailing newline.

        :returns: the header line without a newline.
        :rtype: str
        """
        return '#:GATKReport.v' + self.version + ':' + str(len(self.tables))

    def write(self, filename):
        """
        Write the report to filename.

        :param str filename: The file name to write.
        """
        with open(filename, 'w') as fh:
            fh.write(str(self))

    def __str__(self):
        """
        A string containing the report.

        Each table is separated by two newlines.
        """
        return self.get_headerstring() + '\n' + \
            '\n\n'.join([str(t) for t in self.tables] + [''])

    def __eq__(self, other):
        """
        Test whether two reports are equal. Note that the order of tables
        is significant, and this function will return false if the tables
        are not in the same order.
        """
        if type(other) is type(self):
            if self.version != other.version:
                return False
            else:
                for s, o in zip(self.tables, other.tables):
                    if not s == o:
                        return False
                else:
                    return True
        else:
            return NotImplemented

class GATKTable:
    """
    A class representing a GATK report table.

    The intended method of interacting with the table is using the data
    attribute, which is a standard Pandas dataframe.

    Attributes

        * :attr:`title` - The title of the table
        * :attr:`description` - A description of the table. May be ``''``.
        * :attr:`data` - A Pandas dataframe containing the table data. Accessing this
            attribute is the primary way to interact with the data.
        * :attr:`typemap` - A dict mapping types to a character to be used
            for constructing the format string.
        * :attr:`precisionmap` - A dict mapping column headers to precision
            to be used for constructing the format string.

    Methods

        * :meth:`fromstring` - Initialize a :class:`GATKTable` from a string.
        * :meth:`parse_fmtstring` - Parse a fmtstring into a dictionary of types
        * :meth:`get_fmtstring` - Get a formatstring by inspecting the dataframe.
        * :meth:`get_colfmts` - Get a list of format strings for the data in 
            the dataframe
        * :meth:`get_titlestring` - Get the title string
        * :meth:`get_datastring` - Format the dataframe into a string
        * :meth:`get_nrows` - Get the number of rows in the dataframe.
        * :meth:`get_ncols` - Get the number of columns in the dataframe.
        * :meth:`write` - Write the table to a filehandle.
    """

    def __init__(self, title, description, data):
        """
        Initialize the table.

        :param str title: the table title
        :param str description: a short description of the table. May be ``''``
        :param pandas.DataFrame data: the data contained by the table

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
        self.typemap = {np.dtype(np.int) : 'd', np.dtype(np.float) : 'f', str : 's', np.dtype(np.object) : 's'}
        """
        A dictionary of `{type : str}` to determine the string for each type.

        Note that in the current implementation this only affects
        :meth:`get_fmtstring`.
        """
        self.precisionmap = {'EmpiricalQuality' : '.4', 'EstimatedQReported' : '.4', 'Errors' : '.2'}
        """
        A dictionary of `{colname : str}` to determine the precision string to prepend
        to each formatstring. This only affects :meth:`get_fmtstring`.
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

        :param list(str) header: The column titles as a list
        :param str fmtstring: The format string to parse
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
        frame along with the strings returned by :func:`get_colfmts` to
        create the format string

        which is formatted like ``#:GATKTable:ncol:nrow:%f:%f:%f:;``.

        :return: The format string
        :rtype: str

        """
        fmtlist = ['#', 'GATKTable', str(self.get_ncols()), str(self.get_nrows())]
        for f in self.get_colfmts():
            fmtlist.append(f)
        fmtlist.append(';')
        return ':'.join(fmtlist)

    def get_colfmts(self):
        """
        Get the formats for each column in the dataframe.

        Uses :attr:`typemap` and :attr:`precisionmap` to
        create the strings. This will return something like
        ``['%s', '%s', '%.2f', '%.4f']``

        :return: format strings for each column
        :rtype: list(str)
        """
        fmtlist = []
        unindexed = self.data.reset_index()
        types = unindexed.dtypes
        heads = unindexed.columns.to_list()
        for t,h in zip(types, heads):
            colfmtstr = '%' + self.precisionmap.get(h,'') + self.typemap[t]
            fmtlist.append(colfmtstr)
        return fmtlist


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
        unindexed = self.data.reset_index()
        datawidths = unindexed.apply(lambda x: x.astype('str').str.len().max()).to_numpy()
        headwidths = unindexed.columns.str.len().to_numpy()
        header = unindexed.columns.to_list()
        colwidths = np.maximum(datawidths, headwidths)
        fmtstrings = self.get_colfmts()
        datastr = '  '.join([h.ljust(colwidths[i]) for i,h in enumerate(header)])
        for row in unindexed.itertuples(index = False):
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
        return self.data.reset_index().shape[0]

    def get_ncols(self):
        """
        Get the number of columns in the dataframe.

        :return: The number of columns in the dataframe.
        :rtype: int
        """
        return self.data.reset_index().shape[1]

    def write(self, filehandle):
        """
        Write the table to a filehandle.

        Adds a newline to the end that isn't present in :meth:`__str__`

        :param filehandle: The file object to write the table to.
        :type filehandle: File Object
        :return: Whatever the file object's write method returns.
            Usually the number of characters written.
        :rtype: Usually int

        """
        return filehandle.write(str(self) + '\n')

    def __str__(self):
        return self.get_fmtstring() + '\n' + self.get_titlestring() + '\n' + self.get_datastring()

    def __repr__(self):
        return self.get_fmtstring() + '\n' + self.get_titlestring() + '\n' + repr(self.data)

    def __eq__(self, other):
        if type(other) is type(self):
            return self.title == other.title and \
                self.description == other.description and \
                self.data.equals(other.data)
        else:
            return NotImplemented

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
    2      6                  ReadGroup, EventType, EmpiricalQuality,
                              EstimatedQReported, Observations, Errors
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

    def __init__(self, tables, version = '1.1'):
        """
        Initialize the recalibration report from a list of :class:`GATKTable` s.

        This initializes the tables with some data wrangling to set indices
        and data types properly. If you don't want this, use a 
        :class:`.GATKReport` instead.
        """
        super().__init__(tables, version)

        if len(self.tables) != 5:
            raise ValueError(f"""A RecalibrationReport should have 5 tables.
                This report contains {len(self.tables)}.""")

        #self.data[0] is argument / value
        assert self.tables[0].title == 'Arguments'
        #self.data[1] is quantization map
        assert self.tables[1].title == 'Quantized'
        #self.data[2] is RG
        assert self.tables[2].title == 'RecalTable0'
        #self.data[3] is RG / reportedqual
        assert self.tables[3].title == 'RecalTable1'
        #self.data[4] is RG / reportedqual / covariate (Cycle OR Context)
        assert self.tables[4].title == 'RecalTable2'
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

    def __str__(self):
        """
        GATKReports list CovariateValue before CovariateName, which
        doesn't make much sense for organizing the index. So we switch
        the columns before printing and switch them back after we get
        the string.
        """
        self.tables[4].data = self.tables[4].data.swaplevel('CovariateValue','CovariateName')
        selfstr = super().__str__()
        self.tables[4].data = self.tables[4].data.swaplevel('CovariateValue','CovariateName')
        return selfstr
