from typing import Optional

from hypothesis import strategies
from hypothesis.extra.numpy import arrays
import numpy as np
from hypothesis.searchstrategy import SearchStrategy

from tests.strategies.utils import floats


# TODO: is this a good way to make a strategy?
def finite_nonnegative_floats(max_value: Optional[float] = None
                              ) -> SearchStrategy:
    return strategies.floats(min_value=0.,
                             max_value=max_value,
                             allow_nan=False,
                             allow_infinity=False)


def positive_floats(max_value: Optional[float] = None) -> SearchStrategy:
    return strategies.floats(min_value=0.,
                             max_value=max_value,
                             allow_nan=False,
                             allow_infinity=False).filter(lambda x: x != 0)


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

non_zero_small_floats = small_floats.filter(lambda x: x != 0)
