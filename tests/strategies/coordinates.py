from functools import partial

from hypothesis import strategies
from hypothesis.extra.numpy import arrays
import numpy as np

arrays_sizes = strategies.integers(min_value=1,
                                   max_value=100)

positive_floats = strategies.floats(min_value=0.,
                                    max_value=1e10).filter(lambda x: x > 0.)
positive_floats_w_lower_limit = partial(strategies.floats,
                                        max_value=1e14)

positive_floats_arrays = arrays(dtype=float,
                                shape=arrays_sizes,
                                elements=positive_floats)

positive_floats_arrays_w_size = partial(arrays,
                                        dtype=float,
                                        elements=positive_floats)

triangle_angles = strategies.floats(min_value=0.,
                                    max_value=np.pi).filter(
        lambda x: x != 0. and x != np.pi)

cos_sin_values = strategies.floats(min_value=-1.,
                                   max_value=1.)
