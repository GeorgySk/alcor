import math
import pytest

from tests import strategies
from tests.test_sphere_stars import UNIVERSE_AGE
from tests.utils import example


@pytest.fixture(scope='function')
def halo_birth_initial_time() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def halo_stars_formation_time(halo_birth_initial_time: float) -> float:
    return example(strategies.nonnegative_floats(
            max_value=UNIVERSE_AGE - halo_birth_initial_time))


@pytest.fixture(scope='function')
def bin_initial_time() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def time_increment() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def age() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def formation_rate_parameter() -> float:
    # TODO: is it a good way to filter out zero?
    return example(strategies.floats.filter(lambda x: x != 0))


@pytest.fixture(scope='function')
def max_formation_rate() -> float:
    return example(strategies.nonnegative_floats())


@pytest.fixture(scope='function')
def thick_disk_birth_initial_time() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def exponent() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def min_mass() -> float:
    return example(strategies.nonnegative_floats())


@pytest.fixture(scope='function')
def max_mass(min_mass) -> float:
    return example(strategies.nonnegative_floats()
                   .filter(lambda x: x > min_mass))


@pytest.fixture(scope='function')
def star_formation_rate_param() -> float:
    return example(strategies.floats.filter(lambda x: x != 0))


@pytest.fixture(scope='function')
def thin_disk_age_gyr(star_formation_rate_param) -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE)
                   .filter(
        lambda x: x != 0 and math.exp(-x / star_formation_rate_param) != 1.))


@pytest.fixture(scope='function')
def sigma() -> float:
    return example(strategies.floats)
