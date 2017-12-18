import math
import random
from typing import Callable

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope='function')
def stars_without_theta() -> pd.DataFrame:
    return pd.DataFrame(dict(something=[1, 2, 3]))


@pytest.fixture(scope='function')
def angle_covering_sector() -> float:
    return random.uniform(-2. * math.pi, 2. * math.pi)


@pytest.fixture(scope='function')
def generator() -> Callable[[float, float, float], np.ndarray]:
    return np.random.uniform
