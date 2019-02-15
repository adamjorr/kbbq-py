#!usr/bin/env python3
import numpy as np
import pandas as pd

class GATKReport:
    """
    A base class representing a GATK Report file that contains many tables.

    Variables: 
        header - a header line (such as #:GATKReport.v1.1:5)
        tables - a list of GATKTable objects

    Methods:
        write - write the report to a file

    """

    def __init__(self, filename):
        """Initialize by reading a file specified by filename."""
        with open(filename) as fh:
            self.header = fh.readline()
            table_strings = fh.read().split('\n\n')
            self.tables = [GATKTable(s) for s in table_strings if s != '']

    def write(self, filename):
        """Write the report to filename."""
        with open(filename, 'w') as fh:
            fh.write(self.header)
            for t in self.tables:
                t.write(fh)
                fh.write('\n')

    def __str__():
        pass #TODO

class GATKTable:
    """
    A class representing a GATK report table.

    The intended method of interacting with the table is using the data
    attribute, which is a standard Pandas dataframe.

    Variables:
        format - the format string representing the data in the table.
            This is like: #:GATKTable:ncol:nrow:%f:%f:%f:;
            It is ':' delimited and ';' terminated.
        ncols - The number of columns in the table
        nrows - The number of rows in the table
        fmtstrings - A list of the format strings for each column in the table
        name - The entire line specifying the name of the table. It is
            formatted like '#:GATKTable:title:subtitle'. subtitle can be
            an empty string.
        title - The title of the table
        subtitle - The subtitle of the table
        header - A list containing the title of each column
        data - A Pandas dataframe containing the table data. Accessing this
            attribute is the primary way to interact with the data.

    Methods:
        write - Write the table to a filehandle.
    """

    def __init__(self, tablestring):
        """Initialize the table from a string"""
        rows = tablestring.splitlines()
        self.format = rows[0]
        splitfmt = self.format.split(':')
        self.ncols = splitfmt[2]
        self.nrows = splitfmt[3]
        self.fmtstrings = splitfmt[4:-1]
        self.name = rows[1]
        self.title = self.name.split(':')[2]
        self.header = rows[2].split()
        assert len(self.header) == len(self.fmtstrings)
        strdata = [s.split() for s in rows[3:]]
        d = dict(zip(self.header, zip(*strdata))) #dictionary {colname : coldata}
        self.data = pd.DataFrame(d)
        typedict = {}
        for i,h in enumerate(self.header):
            f = self.fmtstrings[i]
            if f.endswith('d'):
                type = np.int64
            elif f.endswith('f'):
                type = np.float64
            else:
                type = None
            if type is not None:
                typedict[h] = type
        self.data = self.data.astype(typedict)

    def write(self, filehandle):
        """Write the table to a filehandle."""
        filehandle.writelines([self.format + '\n'] + [self.name + '\n'])
        datawidths = self.data.apply(lambda x: x.str.len().max()).get(self.header).values
        headwidths = np.array([len(s) for s in self.header])
        colwidths = np.maximum(datawidths, headwidths)
        print(*[h.ljust(colwidths[i]) for i,h in enumerate(self.header)], sep = '  ', file = filehandle)
        for row in self.data.itertuples(index = False):
            formatted = [getattr(row,self.header[i]).ljust(colwidths[i]) if f == '%s' else (f % float(getattr(row,self.header[i]))).rjust(colwidths[i]) for i,f in enumerate(self.fmtstrings)]
            print(*formatted, sep = '  ', file = filehandle)

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return self.format + self.name + repr(self.data)

class RecalibrationReport(GATKReport):
    def __init__(self, filename):
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
