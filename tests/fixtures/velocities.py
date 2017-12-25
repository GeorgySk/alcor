import numpy as np
import pytest

from tests.utils import example
from tests import strategies


@pytest.fixture(scope='function')
def x_values() -> np.ndarray:
    return example(strategies.floats_arrays)


@pytest.fixture(scope='function')
def y_values(x_values: np.ndarray) -> np.ndarray:
    return example(strategies.floats_arrays_w_size(shape=x_values.shape))


@pytest.fixture(scope='function')
def angles(x_values: np.ndarray) -> np.ndarray:
    return example(strategies.floats_arrays_w_size(shape=x_values.shape))
