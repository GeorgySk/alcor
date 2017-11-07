from functools import partial
import math

import numpy as np
import pandas as pd


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
                   solar_galactocentric_distance: float,
                   oort_a_const: float,
                   oort_b_const: float) -> None:
    halo_stars_mask = stars['galactic_disk_type'] == 'halo'
    thin_disk_stars_mask = stars['galactic_disk_type'] == 'thin'
    thick_disk_stars_mask = stars['galactic_disk_type'] == 'thick'

    halo_stars = stars[halo_stars_mask]
    thin_disk_stars = stars[thin_disk_stars_mask]
    thick_disk_stars = stars[thick_disk_stars_mask]

    set_halo_stars_velocities(
            halo_stars,
            u_peculiar_solar_velocity=u_peculiar_solar_velocity,
            v_peculiar_solar_velocity=v_peculiar_solar_velocity,
            w_peculiar_solar_velocity=w_peculiar_solar_velocity)
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
# More details at: "Simulating Gaia performances on white dwarfs" by S.Torres
def set_halo_stars_velocities(stars: pd.DataFrame,
                              *,
                              u_peculiar_solar_velocity: float,
                              v_peculiar_solar_velocity: float,
                              w_peculiar_solar_velocity: float,
                              ro_param: float = 10.5,
                              l_param: float = 5.5,
                              sigo_param: float = 80.,
                              sigm_param: float = 145.,
                              lsr_motion: float = 220.,
                              gamma: float = 3.4) -> None:
    # "Radial and tangential dispersions"
    galactocentric_distances = np.sqrt(np.power(stars['r_cylindrical'], 2)
                                       + np.power(stars['z_coordinate'], 2))
    xx_params = (galactocentric_distances - ro_param) / l_param

    squared_sigo_param = sigo_param ** 2
    squared_sigm_param = sigm_param ** 2
    sigr2_params = (
         squared_sigm_param * (0.5 - np.arctan(xx_params / np.pi))
         + squared_sigo_param)
    dsigr2dr_params = (-(1. / l_param) * squared_sigm_param
                       / ((1. + np.power(xx_params, 2)) * np.pi))
    sigt2_params = (0.5 * lsr_motion ** 2 + (1. - gamma / 2.) * sigr2_params
                    + (galactocentric_distances / 2.) * dsigr2dr_params)
    sigr_params = np.sqrt(sigr2_params)
    sigt_params = np.sqrt(sigt2_params)
    sigr_params = lsr_motion / math.sqrt(2.)
    sigt_params = lsr_motion / math.sqrt(2.)

    # "Spherical coordinates"
    stars_count = stars.shape[0]
    r_spherical_velocities = sigr_params * np.random.normal(size=stars_count)
    theta_spherical_velocities = (sigt_params
                                  * np.random.normal(size=stars_count))
    phi_spherical_velocities = (sigt_params
                                * np.random.normal(size=stars_count))

    # "Cartesian velocities"
    deltas = np.pi - stars['galactic_longitude'] - stars['theta_cylindrical']
    sin_deltas = np.sin(deltas)
    cos_deltas = np.cos(deltas)
    x_velocities = (-cos_deltas * r_spherical_velocities
                    + sin_deltas * phi_spherical_velocities)
    y_velocities = (sin_deltas * r_spherical_velocities
                    + cos_deltas * phi_spherical_velocities)
    z_velocities = theta_spherical_velocities
    sin_galactic_longitudes = np.sin(stars['galactic_longitude'])
    cos_galactic_longitudes = np.cos(stars['galactic_longitude'])

    u_velocities = (z_velocities * cos_galactic_longitudes
                    + x_velocities * sin_galactic_longitudes)
    v_velocities = (-z_velocities * sin_galactic_longitudes
                    + x_velocities * cos_galactic_longitudes)
    w_velocities = y_velocities

    v_velocities -= lsr_motion

    stars['u_velocity'] = u_velocities + u_peculiar_solar_velocity
    stars['v_velocity'] = v_velocities + v_peculiar_solar_velocity
    stars['w_velocity'] = w_velocities + w_peculiar_solar_velocity
