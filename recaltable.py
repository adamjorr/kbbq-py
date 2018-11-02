#!usr/bin/env python3
import numpy as np
import pandas as pd

class GATKReport:
    def __init__(self, filename):
        with open(filename) as fh:
            self.header = fh.readline()
            table_strings = fh.read().split('\n\n')
            self.tables = [GATKTable(s) for s in table_strings]

    def write(filename):
        pass #TODO

    def __str__():
        pass #TODO

class GATKTable:
    def __init__(self, tablestring):
        rows = tablestring.splitlines()
        self.format = rows[0]
        self.name = rows[1]
        self.title = self.name.split(':')[1]
        self.header = rows[2].split()
        strdata = [s.split() for s in rows[3:]]
        d = dict(zip(header, zip(*strdata))) #dictionary {colname : coldata}
        self.data = pd.DataFrame(d)
