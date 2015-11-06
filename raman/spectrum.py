"""
Defines Spectrum class and associated helper functions.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

from .util import *

def lorentzian(x, amplitude, center, width):
	"""
	Evaluate a lorentzian function with the given parameters at x.
	"""

	numerator =  (width**2)
	denominator = (x - center)**2 + width**2
	return amplitude*(numerator/denominator)


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
		self.x_array = linspace(0, max(self.frequencies), points)
		self.lorentzian = self.lorentzian_sum(width)

	def __eq__(self, other):
		if type(self) is not type(other):
			return False
		return self.frequencies == other.frequencies \
			and self.intensities == other.intensities

	def __len__(self):
		return len(self.frequencies)

	def lorentzian_sum(self, width):
		"""
		Constructs a sum of lorentzians with the given width about the spectral points.
		:param width: the width of each lorentzian, in x-axis units
		"""

		lorentzians = []

		for frequency, intensity in zip(self.frequencies, self.intensities):
			lorentzians.append(
				[lorentzian(point, intensity, frequency, width) for point in self.x_array])

		return [sum(i) for i in zip(*lorentzians)]

	def plot(self, ax, **kwargs):
		"""
		Plot the lorentzian representation of the spectrum.
		:param ax: a matplotlib axis object on which to plot.
		:param kwargs: keyword arguments to be passed to matplotlib.Axis.plot
		"""

		ax.plot(self.x_array, self.lorentzian, **kwargs)

	@property
	def integral(self):
		"""
		Compute the numeric integral of the lorentzian fit to the spectrum.
		"""

		interval = self.x_array[1] - self.x_array[0]
		return sum([value*interval for value in self.lorentzian])

	@staticmethod
	def from_csv(csv_file, points=NUMBER_OF_POINTS):
		"""
		Create a spectrum from a .csv file of frequency-intensity pairs.
		:param csv_file: the path to a .csv file or an opened file object.
		:param points: the number of points to generate in the resulting spectrum.
		"""

		if isinstance(csv_file, str):
			if not csv_file.endswith('.csv'):
				raise ValueError('Filetype must be .csv to create a spectrum.')

			with open(csv_file, 'r') as f:
				return Spectrum.from_csv(f)

		else:
			lines = [line.strip().split(',') for line in csv_file.readlines()]

			frequencies, intensities = [], []
			for line in lines:
				if len(line) > 1 and is_numeric(line[0]) and is_numeric(line[1]):
					frequencies.append(float(line[0]))
					intensities.append(float(line[1]))

			return Spectrum(frequencies, intensities, points)
