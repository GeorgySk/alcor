import pandas as pd
import pytest


@pytest.fixture(scope='function')
def stars_w_progenitor_masses() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[1., 2.]))


@pytest.fixture(scope='function')
def max_mass() -> float:
    return 1.5


@pytest.fixture(scope='function')
def stars_w_galactic_disk_types() -> pd.DataFrame:
    return pd.DataFrame(dict(galactic_disk_type=['thin', 'thick', 'halo']))


@pytest.fixture(scope='function')
def solar_metallicity() -> float:
    return 0.01


@pytest.fixture(scope='function')
def subsolar_metallicity() -> float:
    return 0.001
