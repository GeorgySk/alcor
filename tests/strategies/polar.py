from functools import partial

from hypothesis import strategies
from hypothesis.extra.numpy import arrays
import numpy as np

from .utils import positive_floats

GALACTIC_STRUCTURES = ['thin', 'thick', 'halo']

angles = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-2. * np.pi,
                           max_value=2. * np.pi)

array_sizes = strategies.integers(min_value=1,
                                  max_value=10)

galactic_structures = arrays(
        dtype='<U5',
        shape=array_sizes,
        elements=strategies.sampled_from(GALACTIC_STRUCTURES))

nonnegative_floats = strategies.floats(allow_nan=False,
                                       allow_infinity=False,
                                       min_value=0.,
                                       max_value=1e14)

floats_w_lower_limit = partial(strategies.floats,
                               allow_nan=False,
                               allow_infinity=False,
                               max_value=1e14)

positive_floats_arrays = arrays(dtype=float,
                                shape=array_sizes,
                                elements=positive_floats(max_value=1e14))
