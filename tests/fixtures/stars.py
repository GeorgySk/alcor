import pytest

from alcor.models import Star
from tests import strategies
from tests.utils import example


@pytest.fixture(scope='function')
def star() -> Star:
    return example(strategies.stars)
