import random
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
def floats_tuple() -> Tuple[float, ...]:
    # TODO: find out how to get sorted tuple
    return 0.2456, 1.545645


@pytest.fixture(scope='function')
def cooling_time() -> float:
    return random.uniform(0., 12.)


@pytest.fixture(scope='function')
def cooling_time_grid() -> np.ndarray:
    return np.linspace(0., 9., 10)


@pytest.fixture(scope='function')
def interest_parameter_grid() -> np.ndarray:
    return np.random.rand(10)


@pytest.fixture(scope='function')
def row_index() -> int:
    return random.choice(range(9))


@pytest.fixture(scope='function')
def grid() -> np.ndarray:
    return np.linspace(0., 1., 100)


@pytest.fixture(scope='function')
def metallicity() -> np.ndarray:
    return random.choice([0.001, 0.01, 0.03, 0.06])


@pytest.fixture(scope='function')
def grid_metallicities() -> List[float]:
    return [0.001, 0.01, 0.03, 0.06]


@pytest.fixture(scope='function')
def db_to_da_fraction() -> float:
    return example(strategies.fraction_float)


@pytest.fixture(scope='function')
def size() -> int:
    return 10


@pytest.fixture(scope='function')
def mass() -> float:
    return random.uniform(0.4, 1.4)


@pytest.fixture(scope='function')
def greater_mass_cooling_time_grid() -> np.ndarray:
    return np.linspace(1., 10., 10)


@pytest.fixture(scope='function')
def lesser_mass_cooling_time_grid() -> np.ndarray:
    return np.linspace(0., 9., 10)


@pytest.fixture(scope='function')
def greater_mass_interest_parameter_grid() -> np.ndarray:
    return np.linspace(-12., 4., 10)


@pytest.fixture(scope='function')
def lesser_mass_interest_parameter_grid() -> np.ndarray:
    return np.linspace(-15., 1., 10)


@pytest.fixture(scope='function')
def min_mass() -> float:
    return 0.7


@pytest.fixture(scope='function')
def max_mass() -> float:
    return 1.3


@pytest.fixture(scope='function')
def min_row_index() -> int:
    return 2


@pytest.fixture(scope='function')
def max_row_index() -> int:
    return 3


@pytest.fixture(scope='function')
def star_series() -> pd.Series:
    return pd.Series(dict(mass=1.,
                          luminosity=0.,
                          metallicity=0.01,
                          cooling_time=4.))


@pytest.fixture(scope='function')
def color_table() -> Dict[int, pd.DataFrame]:
    return {52400: pd.DataFrame(dict(
                        u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                        luminosity=np.array([-2., -1., 0., 1.]))),
            100000: pd.DataFrame(dict(
                        u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                        luminosity=np.array([-3., -1., 0., 2.])))}


@pytest.fixture(scope='function')
def color() -> str:
    return 'u_ubvri_absolute'


@pytest.fixture(scope='function')
def tracks() -> Dict[int, pd.DataFrame]:
    return {52400: pd.DataFrame(dict(
                u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                luminosity=np.array([-2., -1., 0., 1.]),
                cooling_time=np.array([3., 4., 5., 6.]))),
            100000: pd.DataFrame(dict(
                u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                luminosity=np.array([-3., -1., 0., 2.]),
                cooling_time=np.array([2., 4., 6., 8.])))}


@pytest.fixture(scope='function')
def interest_parameter() -> str:
    return 'u_ubvri_absolute'


@pytest.fixture(scope='function')
def metallicity_grid() -> List[float]:
    return [0.001, 0.01, 0.03, 0.06]


@pytest.fixture(scope='function')
def cooling_sequences() -> Dict[int, Dict[int, pd.DataFrame]]:
    return {1: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
                100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))},
            10: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))},
            30: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))},
            60: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))}}


@pytest.fixture(scope='function')
def stars() -> pd.DataFrame:
    return pd.DataFrame(dict(mass=[1., 1.1],
                             luminosity=[0., 1.],
                             metallicity=[0.01, 0.001],
                             cooling_time=[3., 4.]))


@pytest.fixture(scope='function')
def stars_without_luminosity() -> pd.DataFrame:
    return pd.DataFrame(dict(mass=[1., 1.1],
                             metallicity=[0.01, 0.001],
                             cooling_time=[3., 4.]))


@pytest.fixture(scope='function')
def da_cooling_sequences() -> Dict[int, Dict[int, pd.DataFrame]]:
    return {1: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))},
            10: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))},
            30: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))},
            60: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))}}


@pytest.fixture(scope='function')
def db_cooling_sequences() -> Dict[int, Dict[int, pd.DataFrame]]:
    return {1: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))},
            10: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))},
            60: {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4]))),
                 100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.]),
                    effective_temperature=np.array([1e4, 1.5e4, 2e4, 3e4])))}}


@pytest.fixture(scope='function')
def da_color_table() -> Dict[int, pd.DataFrame]:
    return {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    b_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    v_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    r_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    i_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    j_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
            100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    b_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    v_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    r_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    i_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    j_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))}


@pytest.fixture(scope='function')
def db_color_table() -> Dict[int, pd.DataFrame]:
    return {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    b_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    v_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    r_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    i_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    j_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
            100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    b_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    v_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    r_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    i_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    j_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))}


@pytest.fixture(scope='function')
def one_color_table() -> Dict[int, pd.DataFrame]:
    return {52400: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-10., 0., 10., 20.]),
                    luminosity=np.array([-2., -1., 0., 1.]),
                    cooling_time=np.array([3., 4., 5., 6.]))),
            100000: pd.DataFrame(dict(
                    u_ubvri_absolute=np.array([-15., -5., 8., 20.]),
                    luminosity=np.array([-3., -1., 0., 2.]),
                    cooling_time=np.array([2., 4., 6., 8.])))}


@pytest.fixture(scope='function')
def random_grid() -> np.ndarray:
    floats_list = example(strategies.floats_lists)
    sorted_floats_list = sorted(example(strategies.floats_lists))
    return np.array([sorted_floats_list,
                     floats_list])


@pytest.fixture(scope='function')
def random_grid_points(random_grid: np.ndarray) -> Tuple[Tuple[float, float],
                                                         Tuple[float, float]]:
    index = random.randrange(random_grid.shape[1] - 1)

    return ((random_grid[0, index], random_grid[0, index + 1]),
            (random_grid[1, index], random_grid[1, index + 1]))


@pytest.fixture(scope='function')
def min_slope(random_grid: np.ndarray) -> float:
    x_array = random_grid[0]
    y_array = random_grid[1]

    x_amplitude = x_array[-1] - x_array[0]
    y_min_distance = array_min_distance(y_array)

    return y_min_distance / x_amplitude


@pytest.fixture(scope='function')
def min_term(random_grid: np.ndarray) -> float:
    x_array = random_grid[0]
    y_array = random_grid[1]

    y_max_distance = y_array.max() - y_array.min()
    x_min_distance = array_min_distance(x_array)

    return y_array.min() - abs(x_array[1]) * y_max_distance / x_min_distance


@pytest.fixture(scope='function')
def max_slope(random_grid: np.ndarray) -> float:
    x_array = random_grid[0]
    y_array = random_grid[1]

    x_min_distance = array_min_distance(x_array)
    y_amplitude = y_array.max() - y_array.min()

    return y_amplitude / x_min_distance


@pytest.fixture(scope='function')
def max_term(random_grid: np.ndarray) -> float:
    x_array = random_grid[0]
    y_array = random_grid[1]

    y_max_distance = y_array.max() - y_array.min()
    x_min_distance = array_min_distance(x_array)

    return y_array.max() + abs(x_array[1]) * y_max_distance / x_min_distance


def array_min_distance(array: np.ndarray) -> float:
    sorted_array = np.sort(array)
    shifted_array = np.roll(sorted_array, 1)
    return (sorted_array[1:] - shifted_array[1:]).min()
