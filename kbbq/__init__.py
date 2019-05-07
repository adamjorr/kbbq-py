"""
kbbq
====

A package for reference-free base quality score recalibration.
"""

import kbbq.compare_reads
from kbbq.compare_reads import *
import kbbq.recaltable
from kbbq.recaltable import *

__all__ = ['compare_reads', 'recaltable', 'plot', 'benchmark']
__version__ = '0.0.0'
