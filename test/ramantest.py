"""
Unit tests for raman module. 

Copyright Sean McGrath 2015

Issued under the MIT licenese.
"""

import unittest
import raman

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


if __name__ == '__main__':
    unittest.main()
