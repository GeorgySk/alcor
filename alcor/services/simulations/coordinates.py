from typing import Union

import numpy as np
import pandas as pd

PC_PER_KPC = 1e3


def stars_w_coordinates(*,
                        stars: pd.DataFrame,
                        solar_galactocentric_distance: float,
                        ngp_declination: float = 0.478,
                        theta: float = 2.147,
                        ngp_right_ascension: float = 3.35,
                        kappa: float = 4.74) -> pd.DataFrame:
    """
    More info on conversions at:
    https://physics.stackexchange.com/questions/88663/converting-between-galactic-and-ecliptic-coordinates
    :param stars: stars without coordinates
    :param solar_galactocentric_distance: distance from Sun to Galactic Center
    :param ngp_declination: declination of the North Galactic Pole
    :param theta: see link above, same as BK = 122.9º there
    :param ngp_right_ascension: right ascension of the North Galactic Pole
    :param kappa: km/sec in one au/year
    :return: stars with assigned coordinates and proper motions
    """
    distance_plane_projections = (opposite_triangle_side(
            solar_galactocentric_distance,
            stars['r_cylindrical'].values,
            stars['theta_cylindrical'].values))
    stars['distance'] = np.sqrt(distance_plane_projections ** 2
                                + stars['z_coordinate'] ** 2)
    stars['galactic_longitude'] = get_galactic_longitudes(
            solar_galactocentric_distance=solar_galactocentric_distance,
            r_cylindrical=stars['r_cylindrical'].values,
            thetas_cylindrical=stars['theta_cylindrical'].values,
            distance_plane_projections=distance_plane_projections)
    stars['galactic_latitude'] = get_galactic_latitudes(
            z_coordinates=stars['z_coordinate'].values,
            distance_plane_projections=distance_plane_projections)

    sin_longitude = np.sin(stars['galactic_longitude'])
    cos_longitude = np.cos(stars['galactic_longitude'])
    sin_latitude = np.sin(stars['galactic_latitude'])
    cos_latitude = np.cos(stars['galactic_latitude'])

    velocities_by_proper_motion = 1. / (kappa * stars['distance'] * PC_PER_KPC)

    # TODO: find out if we need to divide by `cos_latitude`
    stars['proper_motion_in_longitude'] = (
        (velocities_by_proper_motion / cos_latitude)
        * (- stars['u_velocity'] * sin_longitude
           + stars['v_velocity'] * cos_longitude))
    stars['proper_motion_in_latitude'] = velocities_by_proper_motion * (
        - stars['u_velocity'] * cos_longitude * sin_latitude
        - stars['v_velocity'] * sin_latitude * sin_longitude
        + stars['w_velocity'] * cos_latitude)
    stars['radial_velocity'] = (
        stars['u_velocity'] * cos_latitude * cos_longitude
        + stars['v_velocity'] * cos_latitude * sin_latitude
        + stars['w_velocity'] * sin_latitude)

    # TODO: find out if we need to multiply l-component by cos(b)
    stars['proper_motion'] = np.sqrt(stars['proper_motion_in_longitude'] ** 2
                                     + stars['proper_motion_in_latitude'] ** 2)

    stars['declination'] = (np.arcsin(
            np.sin(ngp_declination) * sin_latitude
            + np.cos(ngp_declination) * cos_latitude
            * np.cos(theta - stars['galactic_longitude'])))

    stars['right_ascension'] = right_ascensions(
            cos_latitude=cos_latitude,
            sin_latitude=sin_latitude,
            theta=theta,
            galactic_longitudes=stars['galactic_longitude'],
            declinations=stars['declination'],
            ngp_declination=ngp_declination,
            ngp_right_ascension=ngp_right_ascension)

    return stars


def right_ascensions(*,
                     cos_latitude: np.ndarray,
                     sin_latitude: np.ndarray,
                     theta: float,
                     galactic_longitudes: np.ndarray,
                     declinations: np.ndarray,
                     ngp_declination: float,
                     ngp_right_ascension: float) -> np.ndarray:
    result = np.empty(galactic_longitudes.size)

    # TODO: give better names
    # These variables are for conversion from galactic to equatorial
    # coordinates. More info at the link above
    xs = (cos_latitude * np.sin(theta - galactic_longitudes)
          / np.cos(declinations))
    xc = ((np.cos(ngp_declination) * sin_latitude
           - np.sin(ngp_declination) * cos_latitude
             * np.cos(theta - galactic_longitudes)) / np.cos(declinations))

    # TODO: find out what is going on here
    xs_ge_zero_xc_ge_zero_mask = ((xs >= 0.) & (xc >= 0.))
    xs_ge_zero_xc_lt_zero_mask = ((xs >= 0.) & (xc < 0.))
    xs_lt_zero_xc_ge_zero_mask = ((xs < 0.) & (xc >= 0.))
    xs_lt_zero_xc_lt_zero_mask = ((xs < 0.) & (xc < 0.))

    result[xs_ge_zero_xc_ge_zero_mask] = np.arcsin(xs) + ngp_right_ascension
    result[xs_ge_zero_xc_lt_zero_mask] = np.arccos(xc) + ngp_right_ascension
    result[xs_lt_zero_xc_ge_zero_mask] = (2. * np.pi + np.arcsin(xs)
                                          + ngp_right_ascension)
    result[xs_lt_zero_xc_lt_zero_mask] = (np.pi - np.arcsin(xs)
                                          + ngp_right_ascension)

    mask = result > 2. * np.pi
    result[mask] -= 2. * np.pi

    return result


def get_galactic_latitudes(*,
                           z_coordinates: np.ndarray,
                           distance_plane_projections: np.ndarray
                           ) -> np.ndarray:
    # TODO: or use arctan2?
    result = np.arctan(np.abs(z_coordinates / distance_plane_projections))
    return np.copysign(result, z_coordinates)


def get_galactic_longitudes(*,
                            solar_galactocentric_distance: float,
                            r_cylindrical: np.ndarray,
                            thetas_cylindrical: np.ndarray,
                            distance_plane_projections: np.ndarray
                            ) -> np.ndarray:
    result = triangle_angle(adjacent=np.array([solar_galactocentric_distance]),
                            other_adjacent=distance_plane_projections,
                            opposite=r_cylindrical)
    # Unfolding from 0-180 to 0-360
    mask = thetas_cylindrical > np.pi
    result[mask] = 2. * np.pi - result[mask]

    return result


def opposite_triangle_side(adjacent: Union[float, np.ndarray],
                           other_adjacent: Union[float, np.ndarray],
                           enclosed_angle: Union[float, np.ndarray]
                           ) -> Union[float, np.ndarray]:
    return np.sqrt(adjacent ** 2 + other_adjacent ** 2
                   - 2. * adjacent * other_adjacent * np.cos(enclosed_angle))


def triangle_angle(adjacent: np.ndarray,
                   other_adjacent: np.ndarray,
                   opposite: np.ndarray) -> np.ndarray:
    arccos_argument = ((adjacent ** 2 + other_adjacent ** 2 - opposite ** 2)
                       / (2. * adjacent * other_adjacent))

    positive_erroneous_values_mask = arccos_argument > 1.
    arccos_argument[positive_erroneous_values_mask] = 1.

    negative_erroneous_values_mask = arccos_argument < -1.
    arccos_argument[negative_erroneous_values_mask] = -1.

    return np.arccos(arccos_argument)
