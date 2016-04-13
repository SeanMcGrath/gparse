# gparse

A simple python library for parsing, analyzing and displaying spectral and structural data created by the computational chemsitry program Gaussian.

## Main Classes

### gparse.Spectrum

Represents a complete vibrational spectrum defined by a frequency/intensity for each vibrational mode. Handles the construction of Lorentzian fits to the data, and provides many convenience methods for combining/plotting spectra. Can be instantiated directly from a `.log` file created by Gaussian.

### gparse.PeakAssigner

Parses and structures information regarding the physical vibrations associated with each vibrational mode, including the displacement vector for every atom in a simulated molecule. Can be instantiated directly from a `.log` file created by Gaussian.

### gparse.DistanceMatrix

A convenient data structure for accessing and manipulating the distance matrix associated with a molecular configuration. Can be instantiated directly from a `.log` file created by Gaussian.

## License

Distributed under the MIT license.

## Author

Written by Sean McGrath at UMass Amherst, 2015-2016.
