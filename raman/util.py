"""
Utility functions for implementation of raman analysis.

Part of package raman.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

def is_numeric(string):
    """
    Test if a string is purely numeric.
    :param string: the string to test.
    """
    try:
        float(string)
        return True
    except:
        return False

def linspace(start, end, number_of_points):
    """
    Generate a list of floats from start to end containing number_of_points elements.
    clone of NumPy function with same name.
    """
    if start >= end:
        raise ValueError('The starting value must be less than the ending value.')
    if number_of_points < 2:
        raise ValueError('The space must contain at least two points.')
    interval = (end-start)/(number_of_points-1)
    return [start + interval*i for i in range(number_of_points)]

