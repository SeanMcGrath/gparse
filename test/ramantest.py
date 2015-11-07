"""
Unit tests for raman module. 

Copyright Sean McGrath 2015

Issued under the MIT licenese.
"""

import unittest
import raman

def triangular_number(n):
    return n * (n + 1) // 2

class TestUtilities(unittest.TestCase):
	"""
	Test the methods in raman.util
	"""

	def test_is_numeric(self):
		true_tests = ['3', '4.2', '0000003', 4.56, 1e5, '0.00000000']
		false_tests = ['a', 'a3', '3a', '0.0.0']
		for numeric, not_numeric in zip(true_tests, false_tests):
			self.assertTrue(raman.util.is_numeric(numeric))
			self.assertFalse(raman.util.is_numeric(not_numeric))

	def test_linspace(self):
		self.assertEqual(len(raman.util.linspace(0,10,100)), 100)
		self.assertRaises(ValueError, raman.util.linspace, 10, 0, 100)
		self.assertRaises(ValueError, raman.util.linspace, 0, 10, 1)


class TestSpectrum(unittest.TestCase):
	"""
	Test the Spectrum class.
	"""

	def setUp(self):
		frequencies = range(100)
		intensities = range(0, 1000, 10)
		self.spectrum = raman.Spectrum(frequencies, intensities)

	def test_lorentzian_sum(self):
		self.assertEqual(
			len(self.spectrum.lorentzian_sum(20)),
			len(self.spectrum.x_array))

	def test_from_csv(self):

		self.assertRaises(ValueError, raman.Spectrum.from_csv, 'ramantest.py')
		
		test_spectrum = raman.Spectrum.from_csv('test_spectrum.csv')
		with open('test_spectrum.csv') as csv_file:
			test_spectrum2 = raman.Spectrum.from_csv(csv_file)

		self.assertEqual(test_spectrum, test_spectrum2)


class TestDistanceMatrix(unittest.TestCase):
	"""
	Tests for DistanceMatrix class.
	"""

	def setUp(self):
		self.test_matrix = raman.DistanceMatrix.from_csv('test_matrix.csv')

	def test_len(self):
		test_cases = [((), 0), 
					([(1,2), (1,2,3), (1,2,3,4)], 4)]
		for matrix, length in test_cases:
			self.assertEqual(len(raman.DistanceMatrix(matrix)), length)

	def test_rshift(self):
		test = raman.DistanceMatrix([(1,2), (3,4)])
		test2 = raman.DistanceMatrix([(2, 3), (4,5)])
		self.assertEqual(test >> 1, test2)

	def test_from_csv(self):
		self.assertRaises(ValueError, raman.DistanceMatrix.from_csv, 'ramantest.py')

		with open('test_matrix.csv') as csv_file:
			test_matrix2 = raman.DistanceMatrix.from_csv(csv_file)
			self.assertEqual(self.test_matrix, test_matrix2)

	def test_flattened(self):
		flat_matrix = self.test_matrix.flattened

		# Test that length is preserved
		correct_length = sum([len(row) for row in self.test_matrix])
		self.assertEqual(len(flat_matrix), correct_length)

		# Test that order is preserved
		highest_index = len(self.test_matrix)
		for i, j in zip(range(1, highest_index), range(highest_index - 1)):
			from_test = self.test_matrix[i][j]
			from_flattened = flat_matrix[triangular_number(i) + j]
			self.assertEqual(from_test, from_flattened)

	def test_rms_deviation(self):

		# Shifting matrix should change RMS by value of shift
		for i in range(1,10):
			shifted_matrix = self.test_matrix >> i
			self.assertEqual(self.test_matrix.rms_deviation(shifted_matrix), i)

if __name__ == '__main__':
    unittest.main()
