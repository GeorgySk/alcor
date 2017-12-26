import numpy as np

from alcor.services.simulations.velocities import (rotate_vectors,
                                                   halo_stars_velocities,
                                                   disk_stars_velocities,
                                                   VelocityVector)
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
                               peculiar_solar_velocity: VelocityVector,
                               lsr_velocity: float,
                               spherical_velocity_component_sigma: float,
                               gaussian_generator: GaussianGeneratorType
                               ) -> None:
    x_velocities, y_velocities, z_velocities = halo_stars_velocities(
            galactic_longitudes=galactic_longitudes,
            thetas_cylindrical=thetas_cylindrical,
            peculiar_solar_velocity=peculiar_solar_velocity,
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


def test_disk_stars_velocities(r_cylindrical: np.ndarray,
                               thetas: np.ndarray,
                               peculiar_solar_velocity: VelocityVector,
                               gaussian_generator: GaussianGeneratorType,
                               solar_galactocentric_distance: float,
                               oort_a_const: float,
                               oort_b_const: float,
                               velocity_dispersion: VelocityVector,
                               ) -> None:
    x_velocities, y_velocities, z_velocities = disk_stars_velocities(
            thetas_cylindrical=thetas,
            peculiar_solar_velocity=peculiar_solar_velocity,
            solar_galactocentric_distance=solar_galactocentric_distance,
            oort_a_const=oort_a_const,
            oort_b_const=oort_b_const,
            generator=gaussian_generator,
            r_cylindrical=r_cylindrical,
            velocity_dispersion=velocity_dispersion)

    assert isinstance(x_velocities, np.ndarray)
    assert isinstance(y_velocities, np.ndarray)
    assert isinstance(z_velocities, np.ndarray)
    assert x_velocities.size == r_cylindrical.size
    assert y_velocities.size == r_cylindrical.size
    assert z_velocities.size == r_cylindrical.size
