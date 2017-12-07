import pytest

from tests import strategies
from tests.utils import example


UNIVERSE_AGE = 14.


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
