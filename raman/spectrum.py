"""
Defines Spectrum class and associated helper functions.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

from .util import linspace, is_numeric


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

        self.frequencies        = frequencies
        self.intensities        = intensities
        self.x_array            = linspace(min(self.frequencies), max(self.frequencies), points)
        self.lorentzian_width   = width

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        return self.frequencies == other.frequencies \
            and self.intensities == other.intensities

    def __len__(self):
        return len(self.frequencies)

    def __sub__(self, other):
        """
        Subtract the lorentzians of two spectra.
        """

        # Need a common x array
        highest_frequency = max(max(self.frequencies), max(other.frequencies))
        lowest_frequency = min(min(self.frequencies), min(other.frequencies))
        x_array = linspace(lowest_frequency, highest_frequency, self.NUMBER_OF_POINTS)

        # Need shallow copy of this spectrum
        self_copy = self.copy()

        self_copy.x_array = x_array
        other.x_array = x_array

        return [self_val - other_val for self_val, other_val \
            in zip(self_copy.lorentzian(), other.lorentzian())]


    def lorentzian(self, width=None):
        """
        Constructs a sum of lorentzians with the given width about the spectral points.
        :param width: the width of each lorentzian, in x-axis units
        """

        if not width:
            width = self.lorentzian_width

        lorentzians = []

        for frequency, intensity in zip(self.frequencies, self.intensities):
            lorentzians.append(
                [lorentzian(point, intensity, frequency, width) \
                    for point in self.x_array])

        return [sum(i) for i in zip(*lorentzians)]

    def plot(self, axis, **kwargs):
        """
        Plot the lorentzian representation of the spectrum.
        :param ax: a matplotlib axis object on which to plot.
        :param kwargs: keyword arguments to be passed to matplotlib.Axis.plot
        """

        axis.plot(self.x_array, self.lorentzian(), **kwargs)

    def copy(self):
        """
        Create a shallow copy of this spectrum.
        """

        return Spectrum(self.frequencies, self.intensities, len(self.x_array), self.lorentzian_width)

    def integral(self, width=None):
        """
        Compute the numeric integral of the lorentzian fit to the spectrum.
        :param width: width of lorentzians in the spectrum.
        """

        if not width:
            width = self.lorentzian_width

        interval = self.x_array[1] - self.x_array[0]
        return sum([value * interval for value in self.lorentzian()])

    @staticmethod
    def from_csv(csv_file, points=NUMBER_OF_POINTS, width=LORENTZIAN_WIDTH):
        """
        Create a spectrum from a .csv file of frequency-intensity pairs.
        :param csv_file: the path to a .csv file or an opened file object.
        :param points: the number of points to generate in the resulting spectrum.
        """

        if isinstance(csv_file, str):
            if not csv_file.endswith('.csv'):
                raise ValueError('Filetype must be .csv to create a Spectrum.')

            with open(csv_file, 'r') as open_file:
                return Spectrum.from_csv(open_file)

        else:
            lines = [line.strip().split(',') for line in csv_file.readlines()]

            frequencies, intensities = [], []
            for line in lines:
                if len(line) > 1 and is_numeric(line[0]) and is_numeric(line[1]):
                    frequencies.append(float(line[0]))
                    intensities.append(float(line[1]))

            return Spectrum(frequencies, intensities, points, width)
