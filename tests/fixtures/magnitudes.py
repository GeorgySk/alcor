import random
from typing import Tuple, List

import numpy as np
import pytest

from tests.utils import example
from tests import strategies


@pytest.fixture(scope='function')
def float_value() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def floats_tuple() -> Tuple[float, float]:
    return (example(strategies.floats),
            example(strategies.floats))


@pytest.fixture(scope='function')
def other_floats_tuple() -> Tuple[float, float]:
    return (example(strategies.floats),
            example(strategies.floats))


@pytest.fixture(scope='function')
def cooling_time() -> float:
    return random.uniform(0., 12.)


@pytest.fixture(scope='function')
def cooling_time_grid() -> np.ndarray:
    return np.random.uniform(low=0.,
                             high=12.,
                             size=10)


@pytest.fixture(scope='function')
def interest_parameter_grid() -> np.ndarray:
    return np.random.rand(10)


@pytest.fixture(scope='function')
def row_index() -> int:
    return random.choice(range(10))


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
