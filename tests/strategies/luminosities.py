from functools import partial
from itertools import repeat
from typing import (Any,
                    Callable)

from hypothesis import strategies
from hypothesis.extra.numpy import arrays
from hypothesis.extra.pandas import (data_frames,
                                     columns,
                                     range_indexes)
from hypothesis.searchstrategy import SearchStrategy
import numpy as np
import pandas as pd

from alcor.models.star import GalacticStructureType

GALACTIC_STRUCTURES = [attribute
                       for attribute in dir(GalacticStructureType)
                       if not attribute.startswith('_')]

ASSIGNED_METALLICITIES = [0.001, 0.01]

dataframes_w_galactic_structure_types = data_frames(
        columns=columns(['galactic_structure_type'],
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
