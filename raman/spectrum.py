"""
Defines Spectrum class and associated helper functions.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

from .util import *

DEFAULT_NUMBER_OF_POINTS = 1000

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

	def __init__(self, frequencies, intensities, points=DEFAULT_NUMBER_OF_POINTS):
		"""
		Constructor.
		:param frequencies: a list of frequencies (floats)
		:param intensities: a list of intensities (floats)
		:param points: The number of points in the spectrum for plotting.
		"""

		self.frequencies = frequencies
		self.intensities = intensities
		self.x_array = linspace(0, max(self.frequencies), points)

	def lorentzian_sum(self, width):
		"""
		Constructs a sum of lorentzians with the given width about the spectral points.
		"""
		lorentzians = []
		for frequency, intensity in zip(self.frequencies, self.intensities):
			l = []
			for point in self.x_array:
				l.append(lorentzian(point, intensity, frequency, width))
			lorentzians.append(l)
		# lorentzians = [[lorentzian(point, amp, freq, width) for freq, amp in zip(self.frequencies, self.intensities)]  
		# 				for point in self.x_array]
		return [sum(i) for i in zip(*lorentzians)]

	@staticmethod
	def from_csv(csv_file, points=DEFAULT_NUMBER_OF_POINTS):
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
