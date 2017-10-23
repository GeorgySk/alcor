import pytest

from tests.utils import example
from tests import strategies


@pytest.fixture(scope='function')
def float_value() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def float_value_str(float_value: float) -> str:
    return str(float_value)


@pytest.fixture(scope='function')
def non_float_string() -> str:
    return example(strategies.non_float_strings)
