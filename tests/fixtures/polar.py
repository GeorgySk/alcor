import math
import random

import pandas as pd
import pytest


@pytest.fixture(scope='function')
def stars_without_theta() -> pd.DataFrame:
    return pd.DataFrame(dict(something=[1, 2, 3]))


@pytest.fixture(scope='function')
def angle_covering_sector() -> float:
    return random.uniform(-math.tau, math.tau)
