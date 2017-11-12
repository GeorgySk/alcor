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
def other_floats_tuple() -> Tuple[float, float]:
    return (example(strategies.floats),
            example(strategies.floats))


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
    return 0.2


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
                             metallicity=[0.01, 0.06],
                             cooling_time=[3., 4.]))


@pytest.fixture(scope='function')
def stars_without_luminosity() -> pd.DataFrame:
    return pd.DataFrame(dict(mass=[1., 1.1],
                             metallicity=[0.01, 0.06],
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
