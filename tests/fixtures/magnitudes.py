from typing import (Dict,
                    Tuple,
                    List)

import numpy as np
import pandas as pd
import pytest

from tests.utils import example
from tests import strategies


@pytest.fixture(scope='function')
def float_value() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def same_values_grid() -> Tuple[np.ndarray, np.ndarray]:
    return example(strategies.same_value_grids(
            size=strategies.grid_lengths,
            elements=strategies.floats))


@pytest.fixture(scope='function')
def grid() -> Tuple[np.ndarray, np.ndarray]:
    return example(strategies.grids(
            size=strategies.grid_lengths,
            elements=strategies.floats))


@pytest.fixture(scope='function')
def min_slope(grid: Tuple[np.ndarray, np.ndarray]) -> float:
    x_array = grid[0]
    y_array = grid[1]

    x_amplitude = x_array[-1] - x_array[0]
    y_min_distance = array_min_distance(y_array)

    return y_min_distance / x_amplitude


@pytest.fixture(scope='function')
def min_term(grid: Tuple[np.ndarray, np.ndarray]) -> float:
    x_array = grid[0]
    y_array = grid[1]

    y_max_distance = y_array.max() - y_array.min()
    x_min_distance = array_min_distance(x_array)

    return (y_array.min() - abs(x_array[0]) * y_max_distance / x_min_distance
            - 1.)


@pytest.fixture(scope='function')
def max_slope(grid: Tuple[np.ndarray, np.ndarray]) -> float:
    x_array = grid[0]
    y_array = grid[1]

    x_min_distance = array_min_distance(x_array)
    y_amplitude = y_array.max() - y_array.min()

    return y_amplitude / x_min_distance


@pytest.fixture(scope='function')
def max_term(grid: Tuple[np.ndarray, np.ndarray]) -> float:
    x_array = grid[0]
    y_array = grid[1]

    y_max_distance = y_array.max() - y_array.min()
    x_min_distance = array_min_distance(x_array)

    return (y_array.max() + abs(x_array[0]) * y_max_distance / x_min_distance
            + 1.)


def array_min_distance(array: np.ndarray) -> float:
    sorted_array = np.sort(array)
    shifted_array = np.roll(sorted_array, 1)
    return (sorted_array[1:] - shifted_array[1:]).min()


@pytest.fixture(scope='function')
def cooling_time() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def grid_and_index() -> Tuple[np.ndarray, np.ndarray, float]:
    return example(strategies.grids_and_indices(
            size=strategies.grid_lengths,
            elements=strategies.floats))


@pytest.fixture(scope='function')
def x_array() ->np.ndarray:
    return example(strategies.sorted_arrays_of_unique_values)


@pytest.fixture(scope='function')
def metallicity() -> np.ndarray:
    return example(strategies.metallicities)


@pytest.fixture(scope='function')
def grid_metallicities() -> List[float]:
    return strategies.VALID_METALLICITIES


@pytest.fixture(scope='function')
def db_to_da_fraction() -> float:
    return example(strategies.fraction_floats)


@pytest.fixture(scope='function')
def size() -> int:
    return example(strategies.nonnegative_integers)


@pytest.fixture(scope='function')
def mass() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def greater_mass_grid() -> Tuple[np.ndarray, np.ndarray]:
    return example(strategies.grids(
            size=strategies.grid_lengths,
            elements=strategies.floats))


@pytest.fixture(scope='function')
def lesser_mass_grid() -> Tuple[np.ndarray, np.ndarray]:
    return example(strategies.grids(
            size=strategies.grid_lengths,
            elements=strategies.floats))


@pytest.fixture(scope='function')
def min_and_max_mass() -> Tuple[float, float]:
    return example(strategies.min_and_max_masses(elements=strategies.floats))


@pytest.fixture(scope='function')
def min_row_index() -> int:
    return example(strategies.nonnegative_integers)


@pytest.fixture(scope='function')
def max_row_index() -> int:
    return example(strategies.nonnegative_integers)


@pytest.fixture(scope='function')
def star_series() -> pd.Series:
    return pd.Series(dict(
        mass=example(strategies.nonnegative_floats),
        luminosity=example(strategies.floats),
        metallicity=example(strategies.metallicities),
        cooling_time=example(strategies.nonnegative_floats)))


@pytest.fixture(scope='function')
def color_table() -> Dict[int, pd.DataFrame]:
    return example(strategies.color_tables)


@pytest.fixture(scope='function')
def color() -> str:
    return example(strategies.colors)


@pytest.fixture(scope='function')
def tracks() -> Dict[int, pd.DataFrame]:
    return example(strategies.cooling_tracks)


@pytest.fixture(scope='function')
def interest_parameter() -> str:
    return example(strategies.interest_parameters)


@pytest.fixture(scope='function')
def stars() -> pd.DataFrame:
    return example(strategies.stars_df)


@pytest.fixture(scope='function')
def da_cooling_sequences() -> Dict[int, Dict[int, pd.DataFrame]]:
    return example(strategies.da_cooling_tracks)


@pytest.fixture(scope='function')
def db_cooling_sequences() -> Dict[int, Dict[int, pd.DataFrame]]:
    return example(strategies.db_cooling_tracks)


@pytest.fixture(scope='function')
def da_color_table() -> Dict[int, pd.DataFrame]:
    return example(strategies.cooling_tracks)


@pytest.fixture(scope='function')
def db_color_table() -> Dict[int, pd.DataFrame]:
    return example(strategies.cooling_tracks)


@pytest.fixture(scope='function')
def one_color_table() -> Dict[int, pd.DataFrame]:
    return example(strategies.cooling_tracks)
