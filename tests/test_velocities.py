import numpy as np

from alcor.services.simulations.velocities import (rotate_vectors,
                                                   halo_stars_velocities)
from alcor.types import GaussianGeneratorType


def test_rotate_vectors(*,
                        x_values: np.ndarray,
                        y_values: np.ndarray,
                        angles: np.ndarray) -> None:
    rotated_x_values, rotated_y_values = rotate_vectors(x_values=x_values,
                                                        y_values=y_values,
                                                        angles=angles)

    assert isinstance(rotated_x_values, np.ndarray)
    assert isinstance(rotated_y_values, np.ndarray)
    assert rotated_x_values.size == x_values.size
    assert rotated_y_values.size == y_values.size


def test_halo_stars_velocities(galactic_longitudes: np.ndarray,
                               thetas_cylindrical: np.ndarray,
                               u_peculiar_solar_velocity: float,
                               v_peculiar_solar_velocity: float,
                               w_peculiar_solar_velocity: float,
                               lsr_velocity: float,
                               spherical_velocity_component_sigma: float,
                               gaussian_generator: GaussianGeneratorType
                               ) -> None:
    x_velocities, y_velocities, z_velocities = halo_stars_velocities(
            galactic_longitudes=galactic_longitudes,
            thetas_cylindrical=thetas_cylindrical,
            u_peculiar_solar_velocity=u_peculiar_solar_velocity,
            v_peculiar_solar_velocity=v_peculiar_solar_velocity,
            w_peculiar_solar_velocity=w_peculiar_solar_velocity,
            lsr_velocity=lsr_velocity,
            spherical_velocity_component_sigma=(
                spherical_velocity_component_sigma),
            generator=gaussian_generator)

    assert isinstance(x_velocities, np.ndarray)
    assert isinstance(y_velocities, np.ndarray)
    assert isinstance(z_velocities, np.ndarray)
    assert x_velocities.size == galactic_longitudes.size
    assert y_velocities.size == galactic_longitudes.size
    assert z_velocities.size == galactic_longitudes.size
