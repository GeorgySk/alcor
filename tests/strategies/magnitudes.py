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


VALID_METALLICITIES = [0.001, 0.01, 0.03, 0.06]
DA_METALLICITIES_BY_THOUSAND = [1, 10, 30, 60]
DB_METALLICITIES_BY_THOUSAND = [1, 10, 60]
COLORS_LIST = ['u_ubvri_absolute',
               'b_ubvri_absolute',
               'v_ubvri_absolute',
               'r_ubvri_absolute',
               'i_ubvri_absolute',
               'j_ubvri_absolute']
INTEREST_PARAMETERS_LIST = COLORS_LIST + ['luminosity',
                                          'effective_temperature']


floats = strategies.floats(allow_nan=False,
                           allow_infinity=False,
                           min_value=-1e15,
                           max_value=1e15)
nonnegative_floats = strategies.floats(allow_nan=False,
                                       allow_infinity=False,
                                       min_value=0,
                                       max_value=1e15)
fraction_floats = strategies.floats(allow_nan=False,
                                    allow_infinity=False,
                                    min_value=0.,
                                    max_value=1.)

nonnegative_integers = strategies.integers(min_value=0)
grid_lengths = strategies.integers(min_value=2,
                                   max_value=20)

sorted_arrays_of_unique_values = arrays(dtype=np.float64,
                                        shape=grid_lengths,
                                        elements=floats,
                                        unique=True)

valid_metallicities = strategies.sampled_from(VALID_METALLICITIES)


@strategies.composite
def grids(draw: Callable[[SearchStrategy], Any],
          size: SearchStrategy,
          elements: SearchStrategy,
          dtype: np.dtype = np.float64) -> Any:
    x_array = draw(arrays(dtype=dtype,
                          shape=size,
                          elements=elements,
                          unique=True).map(np.sort))
    y_array = draw(arrays(dtype=dtype,
                          shape=x_array.size,
                          elements=elements))
    return x_array, y_array


@strategies.composite
def grids_and_indices(draw: Callable[[SearchStrategy], Any],
                      size: SearchStrategy,
                      elements: SearchStrategy,
                      dtype: np.dtype = np.float64) -> Any:
    x_array = draw(arrays(dtype=dtype,
                          shape=size,
                          elements=elements,
                          unique=True).map(np.sort))
    y_array = draw(arrays(dtype=dtype,
                          shape=x_array.size,
                          elements=elements))
    index = draw(strategies.one_of(
                    strategies.just(-2),
                    strategies.integers(min_value=0,
                                        max_value=x_array.size - 2)))
    return x_array, y_array, index


@strategies.composite
def same_value_grids(draw: Callable[[SearchStrategy], Any],
                     size: SearchStrategy,
                     elements: SearchStrategy,
                     dtype: np.dtype = np.float64) -> Any:
    x_array = draw(arrays(dtype=dtype,
                          shape=size,
                          elements=elements,
                          unique=True).map(np.sort))
    y_value = draw(elements)
    y_array = draw(arrays(dtype=dtype,
                          shape=x_array.size,
                          elements=strategies.just(y_value)))
    return x_array, y_array


@strategies.composite
def min_and_max_masses(draw: Callable[[SearchStrategy], Any],
                       elements: SearchStrategy) -> Any:
    min_mass = draw(elements)
    epsilon = np.finfo(np.float64).eps
    max_mass = draw(strategies.floats(min_value=min_mass + epsilon))
    return min_mass, max_mass


def sort_column(df: pd.DataFrame,
                *,
                column: str) -> pd.DataFrame:
    df[column] = df[column].sort_values().values
    return df


sort_luminosity_column = partial(sort_column,
                                 column='luminosity')
sort_cooling_time_column = partial(sort_column,
                                   column='cooling_time')
sort_effective_temperature_column = partial(sort_column,
                                            column='effective_temperature')

color_tables = strategies.dictionaries(
        keys=nonnegative_integers,
        values=data_frames(columns=columns(COLORS_LIST + ['luminosity'],
                                           dtype=float,
                                           unique=True),
                           rows=strategies.tuples(*[floats] * 7),
                           index=range_indexes(min_size=2)
                           ).map(sort_luminosity_column),
        min_size=2)

colors = strategies.sampled_from(COLORS_LIST)

cooling_tracks = strategies.dictionaries(
        keys=nonnegative_integers,
        values=data_frames(columns=columns(INTEREST_PARAMETERS_LIST
                                           + ['cooling_time'],
                                           dtype=float,
                                           unique=True),
                           rows=strategies.tuples(*[floats] * 9),
                           index=range_indexes(min_size=2,
                                               max_size=5)
                           ).map(sort_luminosity_column)
                            .map(sort_cooling_time_column)
                            .map(sort_effective_temperature_column),
        min_size=2)

interest_parameters = strategies.sampled_from(INTEREST_PARAMETERS_LIST)

stars_df = data_frames(
        columns=columns(['mass', 'metallicity', 'cooling_time'],
                        dtype=float),
        index=range_indexes(min_size=2),
        rows=strategies.tuples(nonnegative_floats,
                               strategies.sampled_from(VALID_METALLICITIES),
                               nonnegative_floats))

da_cooling_tracks = strategies.fixed_dictionaries(dict(zip(
        DA_METALLICITIES_BY_THOUSAND, repeat(cooling_tracks))))
db_cooling_tracks = strategies.fixed_dictionaries(dict(zip(
        DB_METALLICITIES_BY_THOUSAND, repeat(cooling_tracks))))
