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
		for string in true_tests:
			self.assertTrue(raman.util.is_numeric(string))

if __name__ == '__main__':
    unittest.main()
