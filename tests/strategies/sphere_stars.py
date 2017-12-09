from typing import Optional

from hypothesis import strategies
from hypothesis.extra.numpy import arrays
import numpy as np

from tests.strategies.utils import floats


# TODO: is this a good way to make a strategy?
def nonnegative_floats(max_value: Optional[float] = None):
    return strategies.floats(min_value=0.,
                             max_value=max_value,
                             allow_nan=False,
                             allow_infinity=False)


arrays_lengths = strategies.integers(min_value=1,
                                     max_value=20)

numpy_arrays = arrays(dtype=np.float64,
                      shape=arrays_lengths,
                      elements=floats)

positive_integers = strategies.integers(min_value=1)
small_nonnegative_integers = strategies.integers(min_value=0,
                                                 max_value=100)

small_floats = strategies.floats(min_value=-1e14,
                                 max_value=1e14,
                                 allow_nan=False,
                                 allow_infinity=False)

fractions = strategies.floats(min_value=0.,
                              max_value=1.,
                              allow_nan=False,
                              allow_infinity=False)
