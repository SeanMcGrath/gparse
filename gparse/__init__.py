"""
Package for parsing, manipulating, analyzing and plotting
raman spectral data. Also provides utilities for correlating
raman data with the molecular structures which produce them.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""


from .spectrum import Spectrum, PeakAssigner, PeakReporter
from .matrix import DistanceMatrix
from .configuration import Configuration

__title__ = 'raman'
__version__ = '0.0.1'
__author__ = 'Sean McGrath'
__license__ = 'MIT'
__copyright__ = 'Copyright 2015 Sean McGrath'
