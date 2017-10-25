import pandas as pd
import pytest


@pytest.fixture(scope='function')
def stars() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[1., 2.]))


@pytest.fixture(scope='function')
def max_mass() -> float:
    return 1.5
