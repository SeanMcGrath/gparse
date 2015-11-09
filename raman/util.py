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
    :param start: starting point of list.
    :param end: ending point of list.
    :param number_of_points: number of points in returned list.
    """
    if start >= end:
        raise ValueError('The starting value must be less than the ending value.')
    if number_of_points < 2:
        raise ValueError('The space must contain at least two points.')
    interval = (end-start)/(number_of_points-1)
    return [start + interval*i for i in range(number_of_points)]

def integrate(x_array, y_array):
    """
    Calculate the numeric integral of a 2D data set via the midpoint rule.
    :param x_array: x data along which to integrate.
    :param y_array: y data to integrate.
    """

    assert len(x_array) == len(y_array)

    i = 0
    integral = 0
    while i < len(x_array) - 2:
        average = (y_array[i] + y_array[i+1]) / 2
        interval = x_array[i+1] - x_array[i]
        integral += average * interval
        i += 1

    return integral
