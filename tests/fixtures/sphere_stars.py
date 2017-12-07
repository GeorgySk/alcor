import numpy as np
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
    return example(strategies.small_floats.filter(lambda x: x != 0))


@pytest.fixture(scope='function')
def thin_disk_age_gyr() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE)
                   .filter(lambda x: x != 0))


@pytest.fixture(scope='function')
def sigma() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def times() -> np.ndarray:
    return example(strategies.numpy_arrays)


@pytest.fixture(scope='function')
def burst_init_time() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def birth_rate() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def burst_birth_rate() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def max_age() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def time_bins_count() -> float:
    return example(strategies.positive_integers)


@pytest.fixture(scope='function')
def burst_age() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def initial_mass_function_parameter() -> float:
    return example(strategies.small_floats)


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
def thin_disk_stars_fraction() -> float:
    return example(strategies.fractions.filter(lambda x: x != 0))


@pytest.fixture(scope='function')
def thick_disk_stars_fraction() -> float:
    return example(strategies.fractions)


@pytest.fixture(scope='function')
def thin_disk_stars_count() -> float:
    return example(strategies.positive_integers)


@pytest.fixture(scope='function')
def initial_mass_function_parameter() -> float:
    return example(strategies.small_floats)


@pytest.fixture(scope='function')
def thick_disk_age() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def thick_disk_sfr_param() -> float:
    return example(strategies.small_floats.filter(lambda x: x != 0))


@pytest.fixture(scope='function')
def halo_stars_fraction() -> float:
    return example(strategies.fractions)


@pytest.fixture(scope='function')
def halo_age() -> float:
    return example(strategies.nonnegative_floats(max_value=UNIVERSE_AGE))


@pytest.fixture(scope='function')
def halo_stars_formation_time() -> float:
    return example(strategies.small_floats)
