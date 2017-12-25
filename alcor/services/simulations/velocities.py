from functools import partial
from typing import Tuple

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
                   oort_b_const: float,
                   generator: GaussianGeneratorType = np.random.normal
                   ) -> None:
    halo_stars_mask = stars['galactic_disk_type'] == 'halo'
    thin_disk_stars_mask = stars['galactic_disk_type'] == 'thin'
    thick_disk_stars_mask = stars['galactic_disk_type'] == 'thick'

    halo_stars = stars[halo_stars_mask]
    thin_disk_stars = stars[thin_disk_stars_mask]
    thick_disk_stars = stars[thick_disk_stars_mask]

    (halo_stars['u_velocity'],
     halo_stars['v_velocity'],
     halo_stars['w_velocity']) = halo_stars_velocities(
            galactic_longitudes=halo_stars['galactic_longitude'].values,
            thetas_cylindrical=halo_stars['theta_cylindrical'].values,
            u_peculiar_solar_velocity=u_peculiar_solar_velocity,
            v_peculiar_solar_velocity=v_peculiar_solar_velocity,
            w_peculiar_solar_velocity=w_peculiar_solar_velocity,
            lsr_velocity=lsr_velocity,
            spherical_velocity_component_sigma=lsr_velocity / np.sqrt(2.),
            generator=generator)

    stars_velocities = partial(
            disk_stars_velocities,
            u_peculiar_solar_velocity=u_peculiar_solar_velocity,
            v_peculiar_solar_velocity=v_peculiar_solar_velocity,
            w_peculiar_solar_velocity=w_peculiar_solar_velocity,
            solar_galactocentric_distance=solar_galactocentric_distance,
            oort_a_const=oort_a_const,
            oort_b_const=oort_b_const,
            generator=generator)

    (thin_disk_stars['u_velocity'],
     thin_disk_stars['v_velocity'],
     thin_disk_stars['w_velocity']) = stars_velocities(
            r_cylindrical=thin_disk_stars['r_cylindrical'],
            thetas_cylindrical=thin_disk_stars['theta_cylindrical'],
            u_velocity_dispersion=u_velocity_std_thin_disk,
            v_velocity_dispersion=v_velocity_std_thin_disk,
            w_velocity_dispersion=w_velocity_std_thin_disk)
    (thick_disk_stars['u_velocity'],
     thick_disk_stars['v_velocity'],
     thick_disk_stars['w_velocity']) = stars_velocities(
            r_cylindrical=thick_disk_stars['r_cylindrical'],
            thetas_cylindrical=thick_disk_stars['theta_cylindrical'],
            u_velocity_dispersion=u_velocity_std_thick_disk,
            v_velocity_dispersion=v_velocity_std_thick_disk,
            w_velocity_dispersion=w_velocity_std_thick_disk)


def disk_stars_velocities(*,
                          r_cylindrical: np.ndarray,
                          thetas_cylindrical: np.ndarray,
                          u_peculiar_solar_velocity: float,
                          v_peculiar_solar_velocity: float,
                          w_peculiar_solar_velocity: float,
                          solar_galactocentric_distance: float,
                          oort_a_const: float,
                          oort_b_const: float,
                          u_velocity_dispersion: float,
                          v_velocity_dispersion: float,
                          w_velocity_dispersion: float,
                          generator: GaussianGeneratorType
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # TODO: find out what it means
    uops = (u_peculiar_solar_velocity
            + ((3. - (2. * r_cylindrical) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * r_cylindrical
            * np.sin(thetas_cylindrical))
    vops = (v_peculiar_solar_velocity
            + ((3. - (2. * r_cylindrical) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * r_cylindrical
            * np.cos(thetas_cylindrical)
            - (oort_a_const - oort_b_const) * solar_galactocentric_distance)

    stars_count = r_cylindrical.size

    u_velocities = (u_velocity_dispersion * generator(size=stars_count)
                    + uops)
    v_velocities = (v_velocity_dispersion * generator(size=stars_count)
                    + vops - u_velocity_dispersion ** 2 / 120.)
    w_velocities = (w_velocity_dispersion * generator(size=stars_count)
                    + w_peculiar_solar_velocity)

    return u_velocities, v_velocities, w_velocities


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
