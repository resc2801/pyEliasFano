"""
pyEliasFano.

An Elias-Fano structure represents a monotone non-decreasing sequence of n integers from the universe [0 . . . m) occupying 2n+n⌈log2(m/n)⌉ bits.
"""

__version__ = "0.0.8"
__author__ = 'René Schubotz'
__credits__ = 'German Research Centre for Artificial Intelligence'

from .EliasFano import *
from .UniformlyPartitionedEliasFano import *
