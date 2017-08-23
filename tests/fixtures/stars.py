from typing import List

import pytest

from alcor.models import Star
from tests import strategies
from tests.utils import example


@pytest.fixture(scope='function')
def star() -> Star:
    return example(strategies.defined_stars)


@pytest.fixture(scope='function')
def src_stars() -> List[Star]:
    return example(strategies.defined_stars_lists)


@pytest.fixture(scope='function')
def dst_stars() -> List[Star]:
    return example(strategies.defined_stars_lists)
