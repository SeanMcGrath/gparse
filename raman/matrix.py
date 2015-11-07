"""
Defines Matrix class for representation of molecular distance matrices.

Part of raman package.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

from math import sqrt
from .util import *


class DistanceMatrix:
	"""
	Represents a matrix of intermolecular distances that
	occur in a molecule.
	"""

	SUPPORTED_UNITS = ('a', 'angstroms', 'nm', 'nanometers')
	DEFAULT_UNIT = 'angstroms'

	def __init__(self, matrix=[], units=DEFAULT_UNIT):
		"""
		Constructor.
		:param matrix: a 2D array where array[x][y] represents
		distance from atom x to atom y.
		"""

		self.__matrix = matrix
		if units not in self.SUPPORTED_UNITS:
			raise ValueError('Given units must be one of ' + str(SUPPORTED_UNITS))
		self.units = units

	def __getitem__(self, key):
		return self.__matrix[key]

	def __eq__(self, other):
		if not type(other) is type(self):
			return False
		if not self.units[0] == other.units[0]:
			return False
		for a, b in zip(self.__matrix, other.__matrix):
			for c, d in zip(a, b):
				if c != d:
					return False

		return True

	def __len__(self):
		if not self.__matrix:
			return 0
		return max([len(row) for row in self.__matrix])

	def __str__(self):
		return str(self.__matrix)

	def __rshift__(self, other):
		"""
		Shift every value in the matrix by some number.
		e.g. [(0, 1), (2, 3)] >> 1 = [(1,2), (3,4)]
		:param other: number to shift matrix by.
		"""

		new_matrix = []
		for row in self.__matrix:
			new_matrix.append([item + other for item in row])

		return DistanceMatrix(new_matrix, self.units)

	@property
	def flattened(self):
		"""
		Get a copy of all the values in the matrix in flattened (1D) form.
		"""

		flattened_matrix = []
		for row in self.__matrix:
			flattened_matrix += row
		return flattened_matrix
	
	def rms_deviation(self, other):
		"""
		Calculate the root mean square deviation between two matrices.
		:param other: the other matrix to compare to this matrix.
		"""

		if len(self) != len(other):
			raise ValueError('Matrices have incompatible dimensions.')
		if self.units != other.units:
			raise ValueError('Matrices to compare must have the same units.')

		squared_differences = [(a-b)**2 for a,b in zip(self.flattened, other.flattened)]
		return sqrt(sum(squared_differences)/len(squared_differences))

	@staticmethod
	def from_csv(csv_file, units=DEFAULT_UNIT):
		"""
		Create a DistanceMatrix from the data in a .csv file.
		:param csv_file: either the (string) path to a .csv file, or
		an open file object pointing to a .csv file.
		:param units: the units to give the created matrix.
		"""

		if isinstance(csv_file, str):
			if not csv_file.endswith('.csv'):
				raise ValueError('Filetype must be .csv to create a DistanceMatrix.')

			with open(csv_file, 'r') as f:
				return DistanceMatrix.from_csv(f)

		lines = [line.strip().split(',') for line in csv_file.readlines()]
		
		# Strip trailing empty strings if they exist
		try:
			[line.remove('') for line in lines]
		except:
			pass
		
		highest_line_length = 0
		matrix = []
		for line in lines:
			if len(line) == highest_line_length + 1:
				if all(is_numeric(item) for item in line):
					matrix.append([float(item) for item in line])
					highest_line_length += 1
			else:
				raise ValueError('Data lines in input file should be in order of ascending length.')

		return DistanceMatrix(matrix, units)