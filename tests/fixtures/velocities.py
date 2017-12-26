from typing import NamedTuple

import numpy as np
import pytest

from alcor.services.simulations.velocities import VelocityVector
from alcor.types import GaussianGeneratorType
from tests.utils import example
from tests import strategies


@pytest.fixture(scope='function')
def x_values() -> np.ndarray:
    return example(strategies.floats_arrays)


@pytest.fixture(scope='function')
def y_values(x_values: np.ndarray) -> np.ndarray:
    return example(strategies.floats_arrays_w_size(shape=x_values.shape))


@pytest.fixture(scope='function')
def angles(x_values: np.ndarray) -> np.ndarray:
    return example(strategies.floats_arrays_w_size(shape=x_values.shape))


@pytest.fixture(scope='function')
def galactic_longitudes() -> np.ndarray:
    return example(strategies.floats_arrays)


@pytest.fixture(scope='function')
def thetas_cylindrical(galactic_longitudes: np.ndarray) -> np.ndarray:
    shape = galactic_longitudes.shape
    return example(strategies.floats_arrays_w_size(shape=shape))


@pytest.fixture(scope='function')
def peculiar_solar_velocity() -> NamedTuple:
    return VelocityVector(u=example(strategies.floats),
                          v=example(strategies.floats),
                          w=example(strategies.floats))


@pytest.fixture(scope='function')
def lsr_velocity() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def spherical_velocity_component_sigma() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def gaussian_generator() -> GaussianGeneratorType:
    return np.random.normal


@pytest.fixture(scope='function')
def r_cylindrical() -> np.ndarray:
    return example(strategies.floats_arrays)


@pytest.fixture(scope='function')
def thetas(r_cylindrical: np.ndarray) -> np.ndarray:
    shape = r_cylindrical.shape
    return example(strategies.floats_arrays_w_size(shape=shape))


@pytest.fixture(scope='function')
def oort_constant_a() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def oort_constant_b() -> float:
    return example(strategies.floats)


@pytest.fixture(scope='function')
def velocity_dispersion() -> NamedTuple:
    return VelocityVector(u=example(strategies.floats),
                          v=example(strategies.floats),
                          w=example(strategies.floats))


@pytest.fixture(scope='function')
def solar_galactocentric_distance() -> float:
    return example(strategies.floats)
