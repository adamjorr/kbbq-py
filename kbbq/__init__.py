"""
kbbq
====

A package for reference-free base quality score recalibration.
"""

from . import compare_reads
from .compare_reads import *
from . import recaltable
from .recaltable import *

__all__ = ['compare_reads','recaltable']
