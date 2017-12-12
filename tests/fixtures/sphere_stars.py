import random
from functools import partial
from typing import Callable

import numpy as np
import pytest

from alcor.services.simulations.sphere_stars import monte_carlo_generator
from tests import strategies
from tests.test_sphere_stars import UNIVERSE_AGE
from tests.utils import example


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


@pytest.fixture(scope='function')
def generator() -> Callable[[float, float], float]:
    return random.uniform


@pytest.fixture(scope='function')
def initial_mass_function(initial_mass_function_parameter: float
                          ) -> Callable[[float], float]:
    def function(x: float,
                 exponent: float) -> float:
        return x ** exponent

    return partial(function,
                   exponent=initial_mass_function_parameter)


# TODO: find out what to do with this
@pytest.fixture(scope='function')
def initial_mass_generator(generator: Callable[[float, float], float],
                           min_mass: float,
                           max_mass: float,
                           initial_mass_function: Callable[[float], float],
                           initial_mass_function_parameter: float) -> float:
    max_y = (initial_mass_function(min_mass)
             if initial_mass_function_parameter < 0.
             else initial_mass_function(max_mass))
    return partial(monte_carlo_generator,
                   function=initial_mass_function,
                   min_x=min_mass,
                   max_x=max_mass,
                   max_y=max_y,
                   generator=generator)
