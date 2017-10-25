import numpy as np
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


@pytest.fixture(scope='function')
def rightmost_mass() -> float:
    return 5.


@pytest.fixture(scope='function')
def rightmost_time() -> float:
    return 0.1


@pytest.fixture(scope='function')
def masses(rightmost_mass: float) -> np.ndarray:
    return rightmost_mass + 5. * np.random.rand(10)


@pytest.fixture(scope='function')
def metallicities() -> np.ndarray:
    return np.random.choice([0.001, 0.01], size=10)


@pytest.fixture(scope='function')
def stars_without_cooling_times() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[1., 2.],
                             metallicity=[0.001, 0.01],
                             birth_time=[3., 5.]))


@pytest.fixture(scope='function')
def stars_w_cooling_time() -> pd.DataFrame:
    return pd.DataFrame(dict(cooling_time=[-1., 1.]))


@pytest.fixture(scope='function')
def progenitor_masses() -> np.ndarray:
    return np.array([2., 4., 10.])


@pytest.fixture(scope='function')
def stars_without_masses() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[2., 4., 10.]))


@pytest.fixture(scope='function')
def stars_w_masses() -> pd.DataFrame:
    return pd.DataFrame(dict(mass=[1., 2.]))


@pytest.fixture(scope='function')
def main_sequence_stars() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[0.5, 1., 2.],
                             galactic_disk_type=['thin', 'thick', 'halo'],
                             birth_time=[3., 5., 1.]))
