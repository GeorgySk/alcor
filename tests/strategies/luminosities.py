from hypothesis import strategies
from hypothesis.extra.numpy import arrays
from hypothesis.extra.pandas import (data_frames,
                                     columns,
                                     range_indexes)
import numpy as np

from .utils import positive_floats

GALACTIC_STRUCTURES = ['thin', 'thick', 'halo']

dataframes_w_galactic_structure_types = data_frames(
        columns=columns(['galactic_disk_type'],
                        dtype='<U5'),
        rows=strategies.tuples(strategies.sampled_from(GALACTIC_STRUCTURES)),
        index=range_indexes(min_size=1,
                            max_size=10))

small_positive_ints = strategies.integers(min_value=1,
                                          max_value=10)

positive_floats_arrays = arrays(dtype=np.float64,
                                shape=small_positive_ints,
                                elements=positive_floats(max_value=1e14))
