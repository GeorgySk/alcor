from functools import partial

from hypothesis import strategies
from hypothesis.extra.numpy import arrays
from hypothesis.extra.pandas import (data_frames,
                                     columns,
                                     range_indexes)
import numpy as np

from alcor.models.star import GalacticDiskType

GALACTIC_STRUCTURES = [attribute
                       for attribute in dir(GalacticDiskType)
                       if not attribute.startswith('_')]

ASSIGNED_METALLICITIES = [0.001, 0.01]

dataframes_w_galactic_structure_types = data_frames(
        columns=columns(['galactic_disk_type'],
                        dtype='<U5'),
        rows=strategies.tuples(strategies.sampled_from(GALACTIC_STRUCTURES)),
        index=range_indexes(min_size=1,
                            max_size=10))

positive_floats = strategies.floats(min_value=0.,
                                    max_value=1e14,
                                    allow_nan=False,
                                    allow_infinity=False).filter(lambda x:
                                                                 x != 0)

small_positive_ints = strategies.integers(min_value=1,
                                          max_value=10)

positive_floats_arrays = arrays(dtype=np.float64,
                                shape=small_positive_ints,
                                elements=positive_floats)

metallicities = partial(
        arrays,
        dtype=np.float64,
        elements=strategies.sampled_from(ASSIGNED_METALLICITIES))
