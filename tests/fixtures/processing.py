import pytest
from tests import strategies

from tests.utils import example


@pytest.fixture(scope='function')
def filtration_method() -> str:
    return example(strategies.filtration_methods)
