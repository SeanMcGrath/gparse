"""
Utility functions for implementation of raman analysis.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

def is_numeric(string):
    try:
        float(string)
        return True
    except:
        return False

def linspace(start, end, number_of_points):
	"""
	Generate a list of flaots from start to end containing number_of_points elements.
	clone of NumPy function with same name.
	"""

	interval = (end-start)/(number_of_points-1)
	return [start + interval*i for i in range(number_of_points)]
