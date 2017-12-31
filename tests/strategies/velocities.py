from functools import partial

from hypothesis import strategies
from hypothesis.extra.numpy import arrays

arrays_sizes = strategies.integers(min_value=0,
                                   max_value=100)

floats_arrays = arrays(dtype=float,
                       shape=arrays_sizes)

floats_arrays_w_size = partial(arrays,
                               dtype=float)
