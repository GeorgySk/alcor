from hypothesis import strategies
from hypothesis.extra.numpy import arrays
import numpy as np

from .utils import positive_floats

angles = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-2. * np.pi,
                           max_value=2. * np.pi)

array_sizes = strategies.integers(min_value=1,
                                  max_value=10)

positive_floats_arrays = arrays(dtype=float,
                                shape=array_sizes,
                                elements=positive_floats(max_value=1e14))
