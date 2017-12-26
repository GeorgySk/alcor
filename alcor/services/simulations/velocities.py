from collections import namedtuple
from functools import partial
from typing import Tuple

import numpy as np
import pandas as pd

from alcor.models.star import GalacticDiskType
from alcor.types import GaussianGeneratorType

VelocityVector = namedtuple('VelocityVector', ['u', 'v', 'w'])


def set_velocities(stars: pd.DataFrame,
                   *,
                   thin_disk_velocity_std: VelocityVector(u=32.4,
                                                          v=23.,
                                                          w=18.1),
                   thick_disk_velocity_std: VelocityVector(u=50.,
                                                           v=56.,
                                                           w=34.),
                   peculiar_solar_velocity: VelocityVector(u=-11.,
                                                           v=-12.,
                                                           w=-7.),
                   lsr_velocity: float = -220.,
                   solar_galactocentric_distance: float,
                   oort_a_const: float,
                   oort_b_const: float,
                   generator: GaussianGeneratorType = np.random.normal
                   ) -> None:
    halo_stars_mask = stars['galactic_disk_type'] == GalacticDiskType.halo
    thin_disk_stars_mask = stars['galactic_disk_type'] == GalacticDiskType.thin
    thick_disk_stars_mask = (stars['galactic_disk_type']
                             == GalacticDiskType.thick)

    halo_stars = stars[halo_stars_mask]
    thin_disk_stars = stars[thin_disk_stars_mask]
    thick_disk_stars = stars[thick_disk_stars_mask]

    (halo_stars['u_velocity'],
     halo_stars['v_velocity'],
     halo_stars['w_velocity']) = halo_stars_velocities(
            galactic_longitudes=halo_stars['galactic_longitude'].values,
            thetas_cylindrical=halo_stars['theta_cylindrical'].values,
            peculiar_solar_velocity=peculiar_solar_velocity,
            lsr_velocity=lsr_velocity,
            spherical_velocity_component_sigma=lsr_velocity / np.sqrt(2.),
            generator=generator)

    stars_velocities = partial(
            disk_stars_velocities,
            peculiar_solar_velocity=peculiar_solar_velocity,
            solar_galactocentric_distance=solar_galactocentric_distance,
            oort_a_const=oort_a_const,
            oort_b_const=oort_b_const,
            generator=generator)

    (thin_disk_stars['u_velocity'],
     thin_disk_stars['v_velocity'],
     thin_disk_stars['w_velocity']) = stars_velocities(
            r_cylindrical=thin_disk_stars['r_cylindrical'],
            thetas_cylindrical=thin_disk_stars['theta_cylindrical'],
            velocity_dispersion=thin_disk_velocity_std)
    (thick_disk_stars['u_velocity'],
     thick_disk_stars['v_velocity'],
     thick_disk_stars['w_velocity']) = stars_velocities(
            r_cylindrical=thick_disk_stars['r_cylindrical'],
            thetas_cylindrical=thick_disk_stars['theta_cylindrical'],
            velocity_dispersion=thick_disk_velocity_std)


def disk_stars_velocities(*,
                          r_cylindrical: np.ndarray,
                          thetas_cylindrical: np.ndarray,
                          peculiar_solar_velocity: VelocityVector,
                          solar_galactocentric_distance: float,
                          oort_a_const: float,
                          oort_b_const: float,
                          velocity_dispersion: VelocityVector,
                          generator: GaussianGeneratorType
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # TODO: find out what it means
    uops = (peculiar_solar_velocity.u
            + ((3. - (2. * r_cylindrical) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * r_cylindrical
            * np.sin(thetas_cylindrical))
    vops = (peculiar_solar_velocity.v
            + ((3. - (2. * r_cylindrical) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * r_cylindrical
            * np.cos(thetas_cylindrical)
            - (oort_a_const - oort_b_const) * solar_galactocentric_distance)

    stars_count = r_cylindrical.size

    u_velocities = (velocity_dispersion.u * generator(size=stars_count)
                    + uops)
    # TODO: what is this 120.?
    v_velocities = (velocity_dispersion.v * generator(size=stars_count)
                    + vops - velocity_dispersion.u ** 2 / 120.)
    w_velocities = (velocity_dispersion.w * generator(size=stars_count)
                    + peculiar_solar_velocity.w)

    return u_velocities, v_velocities, w_velocities


# TODO: find out what is going on here
def halo_stars_velocities(*,
                          galactic_longitudes: np.ndarray,
                          thetas_cylindrical: np.ndarray,
                          peculiar_solar_velocity: VelocityVector,
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

    u_velocities += peculiar_solar_velocity.u
    v_velocities += peculiar_solar_velocity.v - lsr_velocity
    w_velocities += peculiar_solar_velocity.w

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
