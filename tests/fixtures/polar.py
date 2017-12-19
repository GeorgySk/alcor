from functools import partial
from typing import (Callable,
                    Tuple)

import numpy as np
import pytest

from tests import strategies
from tests.utils import example


@pytest.fixture(scope='function')
def size() -> int:
    return example(strategies.array_sizes)


@pytest.fixture(scope='function')
def angle_covering_sector() -> float:
    return example(strategies.angles)


@pytest.fixture(scope='function')
def generator() -> Callable[[float, float, float], np.ndarray]:
    return np.random.uniform


@pytest.fixture(scope='function')
def unit_range_generator() -> Callable[[Tuple[int, ...]], np.ndarray]:
    return np.random.rand


@pytest.fixture(scope='function')
def galactic_structures() -> np.ndarray:
    return example(strategies.galactic_structures)


@pytest.fixture(scope='function')
def min_sector_radius() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def max_sector_radius(min_sector_radius) -> float:
    return example(strategies.floats_w_lower_limit(
            min_value=min_sector_radius).filter(lambda x:
                                                x != min_sector_radius))


@pytest.fixture(scope='function')
def halo_core_radius() -> float:
    return example(strategies.positive_floats)


@pytest.fixture(scope='function')
def scale_length() -> float:
    return example(strategies.positive_floats)


@pytest.fixture(scope='function')
def radial_distrib_max() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def squared_min_sector_radius() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def squared_radii_difference() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def r_cylindrical() -> float:
    return example(strategies.positive_floats_arrays)


@pytest.fixture(scope='function')
def sector_radius() -> float:
    return example(strategies.nonnegative_floats)


@pytest.fixture(scope='function')
def scale_height() -> float:
    return example(strategies.positive_floats)


@pytest.fixture(scope='function')
def signs_generator() -> float:
    return partial(np.random.choice,
                   [-1, 1])
