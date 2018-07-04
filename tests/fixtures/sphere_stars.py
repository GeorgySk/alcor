import random
from typing import Callable

import numpy as np
import pytest

from tests import strategies
from tests.utils import example

UNIVERSE_AGE = 14.


@pytest.fixture(scope='session')
def universe_age() -> float:
    return UNIVERSE_AGE


@pytest.fixture(scope='function')
def halo_birth_initial_time(universe_age: float) -> float:
    return example(strategies.finite_nonnegative_floats(
            max_value=universe_age))


@pytest.fixture(scope='function')
def halo_stars_formation_time(halo_birth_initial_time: float,
                              universe_age: float) -> float:
    max_formation_time = universe_age - halo_birth_initial_time
    return example(strategies.finite_nonnegative_floats(
            max_value=max_formation_time))


@pytest.fixture(scope='function')
def bin_initial_time(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


@pytest.fixture(scope='function')
def time_increment(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


# TODO: should this be a fixture? It is not used in test_sphere_stars.py
# but it is used in another fixture
@pytest.fixture(scope='function')
def min_mass() -> float:
    return example(strategies.finite_nonnegative_floats())


@pytest.fixture(scope='function')
def max_mass(min_mass: float) -> float:
    return example(strategies.finite_nonnegative_floats()
                   .filter(lambda x: x > min_mass))


@pytest.fixture(scope='function')
def formation_rate_parameter() -> float:
    return example(strategies.non_zero_small_floats)


@pytest.fixture(scope='function')
def disk_age(universe_age: float) -> float:
    return example(strategies.positive_floats(max_value=universe_age))


@pytest.fixture(scope='function')
def sigma() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def times() -> np.ndarray:
    return example(strategies.numpy_arrays)


@pytest.fixture(scope='function')
def birth_rate() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def burst_birth_rate() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def max_age(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


@pytest.fixture(scope='function')
def time_bins_count() -> float:
    return example(strategies.positive_integers)


@pytest.fixture(scope='function')
def burst_age(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


@pytest.fixture(scope='function')
def max_stars_count() -> float:
    return example(strategies.positive_integers)


@pytest.fixture(scope='function')
def sector_radius_kpc() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def burst_formation_factor() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def mass_reduction_factor() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def initial_mass_function_parameter() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def thick_disk_age(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


@pytest.fixture(scope='function')
def formation_rate_exponent() -> float:
    return example(strategies.non_zero_small_floats)


@pytest.fixture(scope='function')
def stars_count() -> float:
    return example(strategies.small_nonnegative_integers)


@pytest.fixture(scope='function')
def birth_initial_time(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


@pytest.fixture(scope='function')
def burst_initial_time(universe_age: float) -> float:
    return example(
            strategies.finite_nonnegative_floats(max_value=universe_age))


@pytest.fixture(scope='session')
def generator() -> Callable[[float, float], float]:
    return random.uniform


@pytest.fixture(scope='function')
def initial_mass_function(initial_mass_function_parameter: float
                          ) -> Callable[[float], float]:
    def function(x: float) -> float:
        return x ** initial_mass_function_parameter

    return function
