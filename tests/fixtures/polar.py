from typing import Callable

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
