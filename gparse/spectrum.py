"""
Defines Spectrum class and associated helper functions.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

import math
import os
import datetime

try:
    import markdown
except:
    markdown = None

from functools import partial
from .util import linspace, is_numeric, integrate, flatten


def lorentzian(x_value, amplitude, center, width):
    """
    Evaluate a lorentzian function with the given parameters at x_value.
    """

    numerator = (width**2)
    denominator = (x_value - center)**2 + width**2
    return amplitude * (numerator / denominator)


class Spectrum:
    """
    Class to represent one vibrational spectrum.

    Consists of a set of frequencies and their corresponding intensities.
    """

    # Default number of points in the list representation of a spectrum.
    NUMBER_OF_POINTS = 5000

    # Default width of lorentzians fitted to the data.
    LORENTZIAN_WIDTH = 3.3

    def __init__(self, frequencies, intensities, width=LORENTZIAN_WIDTH):
        """
        Constructor.
        :param frequencies: a list of frequencies (floats)
        :param intensities: a list of intensities (floats)
        :param points: The number of points in the spectrum for plotting.
        """

        if len(frequencies) != len(intensities):
            raise ValueError(
                'There must be an equal number of frequencies and intensities.')
        if len(frequencies) < 1:
            raise ValueError(
                'There must be at least one frequency-intensity pair.')
        self.frequencies = frequencies
        self.intensities = intensities
        self.lorentzian_width = width

        # Construct the fit function
        lorentzians = []
        for frequency, intensity in zip(self.frequencies, self.intensities):
            lorentzians.append(
                partial(lorentzian,
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
        Subtract two spectra - returns a function describing the
        diffence spectrum.
        """

        return lambda x: self.fit_function(x) - other.fit_function(x)

    def set_width(self, width):
        self.lorentzian_width = width

    @property
    def fit_function(self):
        """
        Access the computed lorentzian fit to the spectrum.
        """
        return self._fit_function

    def x_array(self, x_min=None, x_max=None, points=NUMBER_OF_POINTS):
        """
        Compute an array of values needed for plotting the
        x-axis of a spectrum.
        """
        x_min = x_min or min(self.frequencies)
        x_max = x_max or max(self.frequencies)
        return linspace(min(self.frequencies), max(self.frequencies), points)

    def as_list(self, points=NUMBER_OF_POINTS):
        """
        Constructs a sum of lorentzians about the spectral points,
        and evaluates it at the given number of points.
        """

        return [self.fit_function(x) for x in self.x_array(points)]

    def plot(self, axis, points=NUMBER_OF_POINTS, stems=False, **kwargs):
        """
        Plot the lorentzian representation of the spectrum.
        l:param ax: a matplotlib axis object on which to plot.
        :param kwargs: keyword arguments to be passed to matplotlib.Axis.plot
        """

        axis.plot(self.x_array(points), self.as_list(points), **kwargs)
        if stems:
            axis.stem(self.frequencies, self.intensities, markerfmt=' ')

    def copy(self):
        """
        Create a shallow copy of this spectrum.
        """

        return Spectrum(self.frequencies, self.intensities, self.lorentzian_width)

    @property
    def integral(self, points=NUMBER_OF_POINTS):
        """
        Compute the numeric integral of the lorentzian fit to the spectrum.
        :param points: number of points to integrate over.
        """

        return integrate(self.x_array(points), self.as_list(points))

    @staticmethod
    def from_csv(csv_file, width=LORENTZIAN_WIDTH):
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

            return Spectrum(frequencies, intensities, width)

    @staticmethod
    def from_log_file(filename, type='raman', width=LORENTZIAN_WIDTH):
        """
        Parse a Gaussian .log file and create a Spectrum.
        :param filename: the path to the .log file
        :param type: 'raman' or 'r' for a raman spectrum, 'ir' or 'infrared'
        for an infrared spectrum
        """

        def _parse_line(line):
            """
            Parse the numbers from a line.
            """
            return [float(item) for item in line.split() if is_numeric(item)]

        with open(filename) as open_file:
            lines = open_file.readlines()

        frequencies = flatten([_parse_line(line.strip())
                               for line in lines if 'Frequencies' in line])
        if type in ('r', 'raman'):
            intensities = flatten([_parse_line(line.strip())
                                   for line in lines if 'Raman Activ' in line])
        elif type in ('ir', 'infrared'):
            intensities = flatten([_parse_line(line.strip())
                                   for line in lines if 'IR Inten' in line])
        else:
            raise ValueError("type must be r, ir, raman, or infrared")

        return Spectrum(frequencies, intensities, width)

    @staticmethod
    def average_function(spectra):
        """
        Compute a function that represents an average over the input spectra.
        :param spectra: an iterable of Spectrum objects
        """
        return lambda x: sum([s.fit_function(x) / len(spectra) for s in spectra])


class SpectralPeak:
    """
    Represents one peak in a spectrum with associated information.
    """

    def __init__(
            self,
            number,
            frequency,
            reduced_mass,
            frc_const,
            ir_intensity,
            raman_activity,
            depolar_p,
            depolar_u):

        self.number = number
        self.frequency = frequency
        self.reduced_mass = reduced_mass
        self.frc_const = frc_const
        self.raman_activity = raman_activity
        self.ir_intensity = ir_intensity
        self.depolar_p = depolar_p
        self.depolar_u = depolar_u
        self.atoms = []

    def __repr__(self):
        return 'Peak #{}: {} 1/cm'.format(self.number, self.frequency)

    def assign(self, heavy_only=False):
        return list(filter(lambda x: x.eigen_sum != 0 and (x.element > 1 or heavy_only is False),
                           sorted(self.atoms, key=lambda x: x.eigen_sum)[::-1]))


class Atom:
    """
    Represents one vibrating atom in a molecule.
    """

    elements = {
        '1': 'hydrogen',
        '6': 'carbon',
        '8': 'oxygen'
    }

    def __init__(self, number, element, x, y, z):

        self.number = number
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.eigen_sum = math.sqrt(x**2 + y**2 + x**2)

    @property
    def el_name(self):
        try:
            return self.elements[self.element]
        except:
            return str(self.element)


class PeakAssigner:

    def __init__(self, log_file, heavy_only=False):

        self.peaks = []
        self.index = 0
        self.heavy_only = heavy_only
        self.done = False

        with open(log_file) as f:
            self.lines = f.readlines()

        for i in range(len(self.lines)):
            if 'and normal coordinates:' in self.lines[i]:
                self.index = i + 1
                break

        while not self.done:
            self.make_peaks()

    def __repr__(self):

        out = ''
        for peak in self.peaks:
            out += str(peak) + '\n'
            for assignment in peak.assign(self.heavy_only):
                out += str(assignment) + '\n'
            out += '\n'

        return out

    @staticmethod
    def parse_floats(line):
        floats = []
        for el in line.split():
            try:
                f = float(el)
                floats.append(f)
            except:
                pass

        return floats

    def make_peaks(self):

        peak_props = [[int(el) for el in self.lines[self.index].split()]]
        self.index += 2

        while 'Atom' not in self.lines[self.index]:
            peak_props.append(self.parse_floats(self.lines[self.index]))
            self.index += 1

        prop_zip = list(zip(*peak_props))
        for props in prop_zip:
            self.peaks.append(SpectralPeak(*props))

        self.index += 1
        while not self.lines[self.index].startswith('        '):
            line = self.lines[self.index]
            if line is '\n':
                self.done = True
                break
            split_line = line.split()
            new_atoms = []
            number = int(split_line[0])
            element = int(split_line[1])
            for i in range(len(prop_zip)):
                x = float(split_line[2 + i * 3])
                y = float(split_line[3 + i * 3])
                z = float(split_line[4 + i * 3])
                self.peaks[-(len(prop_zip) - i)
                           ].atoms.append(Atom(number, element, x, y, z))
            self.index += 1


class PeakReporter:

    def __init__(self, log_file, heavy_only=False):

        self.log_file = log_file
        self.spectrum = Spectrum.from_log_file(log_file)
        self.assignments = PeakAssigner(log_file, heavy_only)

    @staticmethod
    def report_setup(path=None):
        if not path:
            path = 'peak_report'
        try:
            os.mkdir(path)
        except OSError:
            pass
        return path

    def report(self, path=None, plt=None):
        report_path = self.report_setup(path)
        currdir = os.getcwd()
        os.chdir(report_path)

        with open('report.md', 'w') as outfile:
            outfile.write('# Peak Assignment Report\n')

            outfile.write('{} peaks found in {} at {}\n\n'.format(
                len(self.spectrum), self.log_file, datetime.datetime.now()))

            for peak in self.assignments.peaks:
                outfile.write('[{}](#{})  \n'.format(
                    peak.frequency, peak.frequency))

            outfile.write('  \n')

            for peak in self.assignments.peaks:
                outfile.write('## Peak #{}: {} 1/cm <a name="{}"></a>  \n'.format(
                    peak.number, peak.frequency, peak.frequency))
                outfile.write('Raman activity: {}  \n'.format(
                    peak.raman_activity))
                if plt:
                    ax = plt.gca()
                    self.spectrum.plot(ax)

                    plt.plot(self.spectrum.x_array(),
                             [lorentzian(
                                 x,
                                 peak.raman_activity,
                                 peak.frequency,
                                 self.spectrum.lorentzian_width) for x in self.spectrum.x_array()],
                             '--')

                    plt.xlim(peak.frequency * 0.8, peak.frequency * 1.2)
                    plt.ylim(0, peak.raman_activity * 1.2)
                    plot_name = '{}.png'.format(peak.frequency)
                    plt.savefig(plot_name)
                    outfile.write('![{}]({})\n\n'.format(plot_name, plot_name))
                    plt.cla()

                outfile.write(
                    '| Atom # | Element | Sum of vibrational eigenvalues |\n')
                outfile.write(
                    '|--------|---------|--------------------------------|\n')
                for atom in peak.assign():
                    outfile.write('| {} | {} | {} |\n'.format(
                        atom.number, atom.el_name, atom.eigen_sum))
                outfile.write('\n')

        if markdown:
            with open('report.md') as md_file:
                text = md_file.read()
                html = markdown.markdown(
                    text, extensions=['markdown.extensions.tables'])
            with open('report.html', 'w',
                      encoding='utf-8', errors='xmlcharrefreplace') as html_file:
                html_file.write(html)

        else:
            print('Install markdown to convert report to html \n\n pip install markdown')

        os.chdir(currdir)
