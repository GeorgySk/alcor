from hypothesis import strategies
import numpy as np


angles = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-2. * np.pi,
                           max_value=2. * np.pi)

array_sizes = strategies.integers(min_value=1,
                                  max_value=10)

