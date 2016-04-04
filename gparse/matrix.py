"""
Defines Matrix class for representation of molecular distance matrices.

Part of raman package.

Copyright Sean McGrath 2015. Issued under the MIT License.
"""

import re
from math import sqrt
from .util import is_numeric, flatten


class DistanceMatrix:
    """
    Represents a matrix of intermolecular distances that
    occur in a molecule.
    """

    SUPPORTED_UNITS = ('a', 'angstroms', 'nm', 'nanometers')
    DEFAULT_UNIT = 'angstroms'

    # Identifies distance matrix entries in Gaussian .log files
    MATRIX_REGEX = re.compile(r'^\s*\d+\s*[A-Z]\s*(\d+\.\d+\s*)+$')

    def __init__(self, matrix, units=DEFAULT_UNIT):
        """
        Constructor.
        :param matrix: a 2D array where array[x][y] represents
        distance from atom x to atom y.
        """

        self.__matrix = matrix
        if units not in self.SUPPORTED_UNITS:
            raise ValueError('Given units must be one of ' +
                             str(self.SUPPORTED_UNITS))
        self.units = units

    def __getitem__(self, key):
        return self.__matrix[key]

    def __setitem__(self, key, value):
        self.__matrix[key] = value

    def __eq__(self, other):
        try:
            if not self.units[0] == other.units[0]:
                return False
            for self_row, other_row in zip(self, other):
                for self_item, other_item in zip(self_row, other_row):
                    if self_item != other_item:
                        return False

            return True
        except (AttributeError, ValueError):
            return False

    def __sub__(self, other):
        new_matrix = []
        for self_row, other_row in zip(self, other):
            new_matrix.append(
                [self_item - other_item
                 for self_item, other_item in zip(self_row, other_row)])
        return DistanceMatrix(new_matrix, self.units)

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

        return flatten(self.__matrix)

    def rms_deviation(self, other, distance_threshold=None):
        """
        Calculate the root mean square deviation between two matrices.
        :param other: the other matrix to compare to this matrix.
        :param distance_threshold: a distance under which an entry in the matrix must
            fall in order to be factored into the RMS.
        """

        if len(self) != len(other):
            raise ValueError('Matrices have incompatible dimensions.')
        if self.units != other.units:
            raise ValueError('Matrices to compare must have the same units.')

        if distance_threshold:
            squared_differences = \
                [(a - b)**2 for a, b in zip(self.flattened, other.flattened)
                 if abs(a - b) < distance_threshold]

        else:
            squared_differences = \
                [(a - b)**2 for a, b in zip(self.flattened, other.flattened)]

        return sqrt(sum(squared_differences) / len(squared_differences))

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
                raise ValueError(
                    'Filetype must be .csv to create a DistanceMatrix.')

            with open(csv_file, 'r') as open_file:
                return DistanceMatrix.from_csv(open_file)

        lines = [line.strip().split(',') for line in csv_file.readlines()]

        # Strip trailing empty strings if they exist
        try:
            [line.remove('') for line in lines]
        except ValueError:
            pass

        highest_line_length = 0
        matrix = []
        for line in lines:
            if len(line) == highest_line_length + 1:
                if all(is_numeric(item) for item in line):
                    matrix.append([float(item) for item in line])
                    highest_line_length += 1
            else:
                raise ValueError(
                    'Data lines in input file should be in order of ascending length.')

        return DistanceMatrix(matrix, units)

    @staticmethod
    def from_log_file(filename):
        """
        Parse a Gaussian .log file and create a DistanceMatrix.
        :param filename: the path to the .log file
        """

        with open(filename) as open_file:
            lines = open_file.readlines()

        # find the distance matrix lines, clean and split them
        split_lines = []
        for line in lines:
            if re.match(DistanceMatrix.MATRIX_REGEX, line.strip()):
                split_lines.append(line.strip().split())
            # Distance matrix is always terminated by a line containing
            # 'stoich'
            if 'stoich' in line.lower():
                break

        # Build the matrix line by line
        number_atoms = max([int(line[0]) for line in split_lines])

        # Need a matrix to insert values in
        temp_matrix = [[None for i in range(number_atoms)]
                       for i in range(number_atoms)]

        # now we iterate through the lines: the first number in each line
        # is the atom number. We use this to index into the matrix.
        old_index = 0
        column_offset = 0
        for line in split_lines:
            current_column = 0
            current_index = int(line[0]) - 1
            if current_index < old_index:
                column_offset += 5
            if current_index is 1:
                column_offset = 0
            while current_column < len(line) - 2:
                distance = line[current_column + 2]
                temp_matrix[current_index][
                    current_column + column_offset] = distance
                current_column += 1
            old_index = current_index

        # Temp matrix has placeholder elements - remove them
        matrix = [[float(item) for item in row if item] for row in temp_matrix]

        return DistanceMatrix(matrix)
