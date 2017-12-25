from functools import partial
import math
from typing import (Callable,
                    Union,
                    Iterable,
                    Tuple)

import numpy as np
import pandas as pd

from alcor.types import GaussianGeneratorType


def set_velocities(stars: pd.DataFrame,
                   *,
                   u_velocity_std_thin_disk: float = 32.4,
                   v_velocity_std_thin_disk: float = 23.,
                   w_velocity_std_thin_disk: float = 18.1,
                   u_velocity_std_thick_disk: float = 50.,
                   v_velocity_std_thick_disk: float = 56.,
                   w_velocity_std_thick_disk: float = 34.,
                   u_peculiar_solar_velocity: float = -11.,
                   v_peculiar_solar_velocity: float = -12.,
                   w_peculiar_solar_velocity: float = -7.,
                   lsr_velocity: float = -220.,
                   solar_galactocentric_distance: float,
                   oort_a_const: float,
                   oort_b_const: float) -> None:
    halo_stars_mask = stars['galactic_disk_type'] == 'halo'
    thin_disk_stars_mask = stars['galactic_disk_type'] == 'thin'
    thick_disk_stars_mask = stars['galactic_disk_type'] == 'thick'

    halo_stars = stars[halo_stars_mask]
    thin_disk_stars = stars[thin_disk_stars_mask]
    thick_disk_stars = stars[thick_disk_stars_mask]

    (halo_stars['u_velocity'],
     halo_stars['v_velocity'],
     halo_stars['w_velocity']) = halo_stars_velocities(
            galactic_longitudes=halo_stars['galactic_longitude'],
            thetas_cylindrical=halo_stars['theta_cylindrical'],
            u_peculiar_solar_velocity=u_peculiar_solar_velocity,
            v_peculiar_solar_velocity=v_peculiar_solar_velocity,
            w_peculiar_solar_velocity=w_peculiar_solar_velocity,
            lsr_velocity=lsr_velocity,
            spherical_velocity_component_sigma=lsr_velocity / np.sqrt(2.),
            generator=np.random.normal)
    set_stars_velocities = partial(
            set_disk_stars_velocities,
            u_peculiar_solar_velocity=u_peculiar_solar_velocity,
            v_peculiar_solar_velocity=v_peculiar_solar_velocity,
            w_peculiar_solar_velocity=w_peculiar_solar_velocity,
            solar_galactocentric_distance=solar_galactocentric_distance,
            oort_a_const=oort_a_const,
            oort_b_const=oort_b_const)
    set_stars_velocities(
            thin_disk_stars,
            u_velocity_dispersion=u_velocity_std_thin_disk,
            v_velocity_dispersion=v_velocity_std_thin_disk,
            w_velocity_dispersion=w_velocity_std_thin_disk)
    set_stars_velocities(
            thick_disk_stars,
            u_velocity_dispersion=u_velocity_std_thick_disk,
            v_velocity_dispersion=v_velocity_std_thick_disk,
            w_velocity_dispersion=w_velocity_std_thick_disk)


def set_disk_stars_velocities(stars: pd.DataFrame,
                              *,
                              u_peculiar_solar_velocity: float,
                              v_peculiar_solar_velocity: float,
                              w_peculiar_solar_velocity: float,
                              solar_galactocentric_distance: float,
                              oort_a_const: float,
                              oort_b_const: float,
                              u_velocity_dispersion: float,
                              v_velocity_dispersion: float,
                              w_velocity_dispersion: float) -> None:
    # TODO: find out what it means
    uops = (u_peculiar_solar_velocity
            + ((3. - (2. * stars['r_cylindrical'])
                     / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * stars['r_cylindrical']
              * np.sin(stars['theta_cylindrical']))
    vops = (v_peculiar_solar_velocity
            + ((3. - (2. * stars['r_cylindrical'])
                     / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * stars['r_cylindrical']
              * np.cos(stars['theta_cylindrical'])
            - (oort_a_const - oort_b_const) * solar_galactocentric_distance)

    stars_count = stars.shape[0]

    stars['u_velocity'] = (u_velocity_dispersion
                           * np.random.normal(size=stars_count) + uops)
    stars['v_velocity'] = (v_velocity_dispersion *
                           np.random.normal(size=stars_count) + vops
                           - u_velocity_dispersion ** 2 / 120.)
    stars['w_velocity'] = (w_velocity_dispersion
                           * np.random.normal(size=stars_count)
                           + w_peculiar_solar_velocity)


# TODO: find out what is going on here
def halo_stars_velocities(*,
                          galactic_longitudes: np.ndarray,
                          thetas_cylindrical: np.ndarray,
                          u_peculiar_solar_velocity: float,
                          v_peculiar_solar_velocity: float,
                          w_peculiar_solar_velocity: float,
                          lsr_velocity: float,
                          spherical_velocity_component_sigma: float,
                          generator: GaussianGeneratorType
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    More details at: "Simulating Gaia performances on white dwarfs" by S.Torres
    """
    stars_count = galactic_longitudes.shape[0]

    r_spherical_velocities = (spherical_velocity_component_sigma
                              * generator(size=stars_count))
    theta_spherical_velocities = (spherical_velocity_component_sigma
                                  * generator(size=stars_count))
    phi_spherical_velocities = (spherical_velocity_component_sigma
                                * generator(size=stars_count))

    deltas = np.pi - galactic_longitudes - thetas_cylindrical

    x_velocities, y_velocities = rotate_vectors(
            x_values=r_spherical_velocities,
            y_values=phi_spherical_velocities,
            angles=deltas)
    x_velocities = -x_velocities  # TODO: why minus?
    z_velocities = theta_spherical_velocities

    v_velocities, u_velocities = rotate_vectors(x_values=x_velocities,
                                                y_values=z_velocities,
                                                angles=galactic_longitudes)
    w_velocities = y_velocities

    u_velocities += u_peculiar_solar_velocity
    v_velocities += v_peculiar_solar_velocity - lsr_velocity
    w_velocities += w_peculiar_solar_velocity

    return u_velocities, v_velocities, w_velocities


def rotate_vectors(*,
                   x_values: np.ndarray,
                   y_values: np.ndarray,
                   angles: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    sin_angles = np.sin(angles)
    cos_angles = np.cos(angles)

    rotated_x_values = cos_angles * x_values - sin_angles * y_values
    rotated_y_values = sin_angles * x_values + cos_angles * y_values

    return rotated_x_values, rotated_y_values
