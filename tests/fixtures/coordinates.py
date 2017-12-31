import pytest

from tests.utils import example
from tests import strategies


@pytest.fixture(scope='function')
def adjacent() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def other_adjacent() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def opposite(adjacent: float,
             other_adjacent: float) -> float:
    min_value = adjacent + other_adjacent
    return example(strategies.positive_floats_w_lower_limit(
            min_value=min_value).filter(lambda x: 1e14 > x > min_value))


@pytest.fixture(scope='function')
def enclosed_angle() -> float:
    return example(strategies.triangle_angles)


@pytest.fixture(scope='function')
def solar_galactocentric_distance() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def single_r_cylindrical() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def theta_cylindrical() -> float:
    return example(strategies.triangle_angles)


@pytest.fixture(scope='function')
def distance_plane_projections() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def z_coordinates() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def cos_latitude() -> float:
    return example(strategies.cos_sin_values)


@pytest.fixture(scope='function')
def sin_latitude() -> float:
    return example(strategies.cos_sin_values)


@pytest.fixture(scope='function')
def theta() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def galactic_longitude() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def declination() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def ngp_declination() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def ngp_right_ascension() -> float:
    return example(strategies.floats)
