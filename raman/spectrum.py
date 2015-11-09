"""
Defines Spectrum class and associated helper functions.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

from .util import linspace, is_numeric, integrate
from functools import partial


def lorentzian(x_value, amplitude, center, width):
    """
    Evaluate a lorentzian function with the given parameters at x_value.
    """

    numerator = (width**2)
    denominator = (x_value - center)**2 + width**2
    return amplitude * (numerator / denominator)


class Spectrum:
    """
    Class to represent one raman spectrum.

    Consists of a set of frequencies and their corresponding intensities.
    """

    NUMBER_OF_POINTS = 1000
    LORENTZIAN_WIDTH = 20

    def __init__(self, frequencies, intensities,
                 points=NUMBER_OF_POINTS, width=LORENTZIAN_WIDTH):
        """
        Constructor.
        :param frequencies: a list of frequencies (floats)
        :param intensities: a list of intensities (floats)
        :param points: The number of points in the spectrum for plotting.
        """

        self.frequencies = frequencies
        self.intensities = intensities
        self.x_array = linspace(min(self.frequencies), max(self.frequencies), points)
        self.lorentzian_width = width

        # Construct the fit function
        lorentzians = []
        for frequency, intensity in zip(self.frequencies, self.intensities):
            lorentzians.append(
                partial(lorentzian, \
                    amplitude=intensity, center=frequency, width=self.lorentzian_width))

        self._fit_function = lambda x: sum([f(x) for f in lorentzians])

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        return self.frequencies == other.frequencies \
            and self.intensities == other.intensities

    def __len__(self):
        return len(self.frequencies)

    def __sub__(self, other):
        """
        Subtract two spectra - returns the difference of the lorentzian fits,
        evaluated at all points in this spectrum's x_array
        """

        difference_function = lambda x: self.fit_function(x) - other.fit_function(x)
        return [difference_function(x) for x in self.x_array]

    @property
    def fit_function(self):
        """
        Access the computed lorentzian fit to the spectrum.
        """
        return self._fit_function

    @property
    def lorentzian(self):
        """
        Constructs a sum of lorentzians with the given width about the spectral points.
        """

        return [self.fit_function(x) for x in self.x_array]

    def plot(self, axis, **kwargs):
        """
        Plot the lorentzian representation of the spectrum.
        :param ax: a matplotlib axis object on which to plot.
        :param kwargs: keyword arguments to be passed to matplotlib.Axis.plot
        """

        axis.plot(self.x_array, self.lorentzian, **kwargs)

    def copy(self):
        """
        Create a shallow copy of this spectrum.
        """

        return Spectrum(self.frequencies, self.intensities, \
            len(self.x_array), self.lorentzian_width)

    @property
    def integral(self):
        """
        Compute the numeric integral of the lorentzian fit to the spectrum.
        :param width: width of lorentzians in the spectrum.
        """

        return integrate(self.x_array, self.lorentzian)

    @staticmethod
    def from_csv(csv_file, points=NUMBER_OF_POINTS, width=LORENTZIAN_WIDTH):
        """
        Create a spectrum from a .csv file of frequency-intensity pairs.
        :param csv_file: the path to a .csv file.
        :param points: the number of points to generate in the resulting spectrum.
        """

        if not csv_file.endswith('.csv'):
            raise ValueError('Filetype must be .csv to create a Spectrum.')

        with open(csv_file, 'r') as open_file:
            lines = [line.strip().split(',') for line in open_file.readlines()]

            frequencies, intensities = [], []
            for line in lines:
                if len(line) > 1 and is_numeric(line[0]) and is_numeric(line[1]):
                    frequencies.append(float(line[0]))
                    intensities.append(float(line[1]))

            return Spectrum(frequencies, intensities, points, width)
