import numpy as np
import pandas as pd
import pytest

from alcor.services.simulations.luminosities import (linear_estimation,
                                                     immutable_array)
from tests import strategies
from tests.utils import example


@pytest.fixture(scope='function')
def stars_w_galactic_structure_types() -> pd.DataFrame:
    return example(strategies.dataframes_w_galactic_structure_types)


@pytest.fixture(scope='session')
def solar_metallicity() -> float:
    return 0.01


@pytest.fixture(scope='session')
def subsolar_metallicity() -> float:
    return 0.001


@pytest.fixture(scope='function')
def rightmost_mass() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def rightmost_time() -> float:
    return example(strategies.positive_floats(max_value=1e14))


@pytest.fixture(scope='function')
def masses(rightmost_mass: float) -> np.ndarray:
    return (np.array([rightmost_mass])
            + example(strategies.positive_floats_arrays))


@pytest.fixture(scope='function')
def metallicities(masses) -> np.ndarray:
    return example(strategies.metallicities(shape=masses.shape))


@pytest.fixture(scope='function')
def stars_without_cooling_times() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[1., 2.],
                             metallicity=[0.001, 0.01],
                             birth_time=[3., 5.]))


@pytest.fixture(scope='function')
def progenitor_masses() -> np.ndarray:
    return example(strategies.positive_floats_arrays)


@pytest.fixture(scope='function')
def stars_without_masses() -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[2., 4., 10.]))


@pytest.fixture(scope='function')
def galactic_disks_types() -> np.ndarray:
    return np.array(['thin', 'thick', 'halo'])


@pytest.fixture(scope='function')
def main_sequence_stars(galactic_disks_types) -> pd.DataFrame:
    return pd.DataFrame(dict(progenitor_mass=[0.5, 1., 2.],
                             galactic_disk_type=galactic_disks_types,
                             birth_time=[3., 5., 1.]))


@pytest.fixture(scope='function')
def model_solar_masses() -> np.ndarray:
    return immutable_array([1.00, 1.50, 1.75, 2.00, 2.25,
                            2.50, 3.00, 3.50, 4.00, 5.00])


@pytest.fixture(scope='function')
def model_solar_times() -> np.ndarray:
    return immutable_array([8.614, 1.968, 1.249, 0.865, 0.632,
                            0.480, 0.302, 0.226, 0.149, 0.088])


@pytest.fixture(scope='function')
def spline(model_solar_masses: np.ndarray,
           model_solar_times: np.ndarray) -> np.ndarray:
    return linear_estimation(x=model_solar_masses,
                             y=model_solar_times)
